/**
 * Implementation of a Rechester-Rosenbluth operator for heat transport.
 *
 * Notes on added physics (relative to the original RR operator):
 *  - We apply a multiplicative correction to the quasi-linear heat diffusivity (original version):
 *      chi_QL -> eta * chi_QL0
 *    where eta is taken from Xiao NF (2025), Eq.(24) (collision/orbit correction).
 *  - We use estimation in Xiao PoP (2024):
 *      D_FL ~ <b^2> L0,
 *      L_S = R q / s,
 *      k_h ~ m/a,  l_cr ~ 1/k_h,
 *      L_K ~ L_S * (k_h^2 D_FL L_S/3)^(-1/3).
 *  - We ignore the Jacobian contribution of eta via current (i.e. via q,s from j_tot).
 *  In SetPartialDiffusionTerm(), L_K, Omega_ce and l_cr are treated constant wrt n,T
 */

#include <cmath>
#include <vector>
#include <limits>

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/HeatTransportRechesterRosenbluth.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"

using namespace DREAM;

// ===================== SOFT-MIN HELPERS (NEW) =====================
// Smooth approximation of min() in chi-space:
//   chiEff = (sum_i chi_i^{-p})^{-1/p}, with p>0 (only positive chi_i contribute).
// This provides a differentiable "hard switch" between competing transport channels.
namespace {
    inline real_t softmin_pos(const std::vector<real_t>& x, const real_t p){
        real_t sum = 0;
        for (real_t xi : x){
            if (xi > 0 && std::isfinite(xi))
                sum += std::pow(xi, -p);
        }
        if (!(sum > 0))
            return 0;
        return std::pow(sum, -1.0/p);
    }

    // For chiEff = (Σ chi^{-p})^{-1/p}, the partial derivative is:
    //   d chiEff / d chi_i = chiEff^(p+1) * chi_i^(-p-1)
    inline real_t dsoftmin_dxi(const real_t xeff, const real_t xi, const real_t p){
        if (!(xeff > 0) || !std::isfinite(xeff)) return 0;
        if (!(xi   > 0) || !std::isfinite(xi))   return 0;
        return std::pow(xeff, p+1.0) * std::pow(xi, -p-1.0);
    }
}
// ==================================================================


/**
 * Constructor.
 */
HeatTransportRechesterRosenbluth::HeatTransportRechesterRosenbluth(
    FVM::Grid *grid,
    FVM::Interpolator1D *dB_B,
    FVM::UnknownQuantityHandler *unknowns,
    RunawayFluid *REFluid,
    real_t L0, real_t Lcr
) : FVM::DiffusionTerm(grid),
    deltaBOverB(dB_B),
    unknowns(unknowns),
    REFluid(REFluid),
    L0(L0), Lcr(Lcr) {

    SetName("HeatTransportRechesterRosenbluth");

    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_j_tot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);

    AddUnknownForJacobian(unknowns, this->id_n_cold);
    AddUnknownForJacobian(unknowns, this->id_T_cold);

    // Soft-min hardness parameter p (NEW): larger p -> closer to hard min.
    this->pSoft = 4.0;

    AllocateDiffCoeff();
}

/**
 * Destructor.
 */
HeatTransportRechesterRosenbluth::~HeatTransportRechesterRosenbluth() {
    if (this->deltaBOverB != nullptr)
        delete this->deltaBOverB;

    delete [] this->dD;

    // (NEW) caches used for Jacobian evaluation
    delete [] this->dBB_cache;
    delete [] this->chiQL_cache;
    delete [] this->chi2_cache;
    delete [] this->chi3_cache;
    delete [] this->chiEff_cache;

    // (NEW) caches for eta-model (NF2025/PoP2024)
    delete [] this->eta_cache;
    delete [] this->LK_cache;
    delete [] this->LKdelta_cache;
}


/**
 * Allocate memory for the differentiation coefficient.
 */
void HeatTransportRechesterRosenbluth::AllocateDiffCoeff() {
    const len_t nr = this->grid->GetNr();
    this->dD = new real_t[nr+1];

    // (NEW) cache dB/B at faces for Jacobian-related evaluations
    this->dBB_cache = new real_t[nr+1];

    // (NEW) caches of chi-branches (faces) for the soft-min Jacobian
    this->chiQL_cache  = new real_t[nr+1];
    this->chi2_cache   = new real_t[nr+1];
    this->chi3_cache   = new real_t[nr+1];
    this->chiEff_cache = new real_t[nr+1];

    // (NEW) caches for eta, LK and LKdelta (faces)
    this->eta_cache     = new real_t[nr+1];
    this->LK_cache      = new real_t[nr+1];
    this->LKdelta_cache = new real_t[nr+1];

    for (len_t ir=0; ir<nr+1; ir++){
        this->dBB_cache[ir] = 0;

        this->chiQL_cache[ir]  = 0;
        this->chi2_cache[ir]   = 0;
        this->chi3_cache[ir]   = 0;
        this->chiEff_cache[ir] = 0;

        this->eta_cache[ir]     = 1.0;
        this->LK_cache[ir]      = 0.0;
        this->LKdelta_cache[ir] = 0.0;
    }
}

/**
 * Called whenever the grid is rebuilt.
 */
bool HeatTransportRechesterRosenbluth::GridRebuilt() {
    this->FVM::DiffusionTerm::GridRebuilt();

    delete [] this->dD;

    // (NEW) delete caches and re-allocate on the new grid
    delete [] this->dBB_cache;
    delete [] this->chiQL_cache;
    delete [] this->chi2_cache;
    delete [] this->chi3_cache;
    delete [] this->chiEff_cache;

    delete [] this->eta_cache;
    delete [] this->LK_cache;
    delete [] this->LKdelta_cache;

    AllocateDiffCoeff();

    return true;
}

/**
 * Rebuild the coefficients for this equation term.
 */
void HeatTransportRechesterRosenbluth::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const real_t *dB_B = this->EvaluateDeltaBOverB(t);
    const len_t nr = this->grid->GetNr();
    const real_t mc2 = Constants::mc2inEV;

    const real_t *ncold = unknowns->GetUnknownData(this->id_n_cold);
    const real_t *Tcold = unknowns->GetUnknownData(this->id_T_cold);
    const real_t *j_tot = unknowns->GetUnknownData(this->id_j_tot);

    FVM::RadialGrid *rg = this->grid->GetRadialGrid();
    const real_t R0 = rg->GetR0();

    // Collision time (thermal) from RunawayFluid (cell centers); we face-interpolate.
    const real_t *tauTh = nullptr;
    if (this->REFluid != nullptr)
        tauTh = this->REFluid->GetElectronCollisionTimeThermal();

    const real_t PREFAC = 3.0 * std::sqrt(2.0*M_PI) * Constants::c * Constants::ec;

    // ===================== q(r) and magnetic shear s(r) (NEW) =====================
    // We estimate s(r) = (r/q) dq/dr on cell centers using q computed from j_tot.
    // This is used in PoP(2024) estimate of the shear length:
    //   L_S = R*q/s
    // Note: we do not include derivatives wrt j_tot in the Jacobian (handled later if desired).
    std::vector<real_t> r_c(nr, 0), q_c(nr, 0), s_c(nr, 0);

    for (len_t i = 0; i < nr; i++) {
        const real_t rf0 = rg->GetR_f(i);
        const real_t rf1 = rg->GetR_f(i+1);
        r_c[i] = 0.5*(rf0 + rf1);
    }

    if (j_tot != nullptr && std::isfinite(R0) && R0 != 0.0) {
        for (len_t i = 0; i < nr; i++) {
            const real_t mu0Ip =
                Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(i, rg, j_tot);
            const real_t qR0_i = rg->SafetyFactorNormalized(i, mu0Ip);
            q_c[i] = qR0_i / R0;
        }

        if (nr >= 2) {
            // one-sided at i=0
            {
                const real_t dr = r_c[1]-r_c[0];
                if (dr != 0 && std::isfinite(dr) && q_c[0] != 0 && std::isfinite(q_c[0])) {
                    const real_t dqdr = (q_c[1]-q_c[0]) / dr;
                    s_c[0] = (r_c[0]/q_c[0]) * dqdr;
                }
            }
            // central in interior
            for (len_t i = 1; i < nr-1; i++) {
                const real_t dr = r_c[i+1]-r_c[i-1];
                if (dr != 0 && std::isfinite(dr) && q_c[i] != 0 && std::isfinite(q_c[i])) {
                    const real_t dqdr = (q_c[i+1]-q_c[i-1]) / dr;
                    s_c[i] = (r_c[i]/q_c[i]) * dqdr;
                }
            }
            // one-sided at i=nr-1
            {
                const real_t dr = r_c[nr-1]-r_c[nr-2];
                if (dr != 0 && std::isfinite(dr) && q_c[nr-1] != 0 && std::isfinite(q_c[nr-1])) {
                    const real_t dqdr = (q_c[nr-1]-q_c[nr-2]) / dr;
                    s_c[nr-1] = (r_c[nr-1]/q_c[nr-1]) * dqdr;
                }
            }
        }
    }
    // ============================================================================

    for (len_t ir = 0; ir < nr+1; ir++) {

        // Face interpolation of T and n
        real_t T=0, n=0;
        if (ir < nr) {
            T += deltaRadialFlux[ir] * Tcold[ir];
            n += deltaRadialFlux[ir] * ncold[ir];
        }
        if (ir > 0) {
            T += (1.0 - deltaRadialFlux[ir]) * Tcold[ir-1];
            n += (1.0 - deltaRadialFlux[ir]) * ncold[ir-1];
        }

        const real_t Theta = T / mc2;

        const real_t B_FSA  = rg->GetFSA_B_f(ir);
        const real_t xiT0   = rg->GetXi0TrappedBoundary_fr(ir);

        // qR0 at this face (use outermost cell for the last face)
        real_t qR0 = 1.0;
        if (j_tot != nullptr) {
            len_t ir_q = ir;
            if (ir_q >= nr && nr > 0) ir_q = nr-1;

            const real_t mu0Ip =
                Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir_q, rg, j_tot);
            qR0 = rg->SafetyFactorNormalized(ir_q, mu0Ip);
        }

        // "D" here is the RR prefactor used in the original implementation (stored for Jacobian).
        const real_t dBB = dB_B[ir];
        const real_t D = PREFAC * dBB*dBB * B_FSA * (1.0 - xiT0*xiT0);
        this->dD[ir] = D;

        // (NEW) store dB/B at faces for possible future Jacobian extensions
        this->dBB_cache[ir] = dBB;

        // Thermal speed vTe and collision frequency nuc (both used in eta and chi3)
        real_t vTe = 0;
        if (Theta > 0 && std::isfinite(Theta))
            vTe = Constants::c * std::sqrt(2.0*Theta);

        real_t nuc = 0;
        if (tauTh != nullptr) {
            real_t tau_face = 0;
            real_t wsum = 0;

            if (ir < nr && std::isfinite(tauTh[ir])) {
                const real_t w = deltaRadialFlux[ir];
                tau_face += w * tauTh[ir];
                wsum += w;
            }
            if (ir > 0 && std::isfinite(tauTh[ir-1])) {
                const real_t w = 1.0 - deltaRadialFlux[ir];
                tau_face += w * tauTh[ir-1];
                wsum += w;
            }

            if (wsum > 0 && tau_face > 0 && std::isfinite(tau_face))
                nuc = 1.0 / tau_face;
        }

        // ===================== eta-model (NEW) =====================
        // Xiao NF (2025), Eq.(24):
        //   eta = 1 / (1 + LKdelta/lambda_mfp + l_cr/rho_e)
        // where
        //   lambda_mfp = vTe / nuc,
        //   rho_e      = vTe / Omega_ce,
        //   LKdelta    = LK * ln[ (sqrt(chi_parallel/chi_perp))*(l_cr/LK) ].
        //
        // Using PoP (2024) classical transport scalings:
        //   chi_parallel ~ v_th^2 / nu,   chi_perp ~ rho_e^2 * nu
        // => sqrt(chi_parallel/chi_perp) = Omega_ce / nu
        // hence
        //   LKdelta = LK * ln[ (Omega_ce/nuc) * (l_cr/LK) ].
        //
        // PoP (2024) QL estimate of LK:
        //   D_FL ~ <b^2> L0,
        //   L_S = R q / s,
        //   k_h ~ m/a,  l_cr ~ 1/k_h,
        //   LK ~ L_S * (k_h^2 D_FL L_S / 3)^(-1/3).
        real_t eta = 1.0;
        real_t LK = 0.0;
        real_t LKdelta = 0.0;

        if (Theta > 0 && vTe > 0 && nuc > 0 && std::isfinite(vTe) && std::isfinite(nuc)
            && j_tot != nullptr && std::isfinite(R0) && R0 != 0.0
            && dBB > 0 && std::isfinite(dBB)) {

            // k_h = m/a,  l_cr = 1/k_h (PoP 2024). Use m=30 as a representative value.
            const real_t a_minor = rg->GetR_f(nr);
            const real_t m_pol = 30.0;
            const real_t kh = (a_minor > 0 && std::isfinite(a_minor)) ? (m_pol / a_minor) : 0.0;
            const real_t lcr = (kh > 0 && std::isfinite(kh)) ? (1.0 / kh) : 0.0;

            // D_FL ~ <b^2> L0 (PoP 2024), with <b^2> ~ (dB/B)^2
            const real_t DFL = dBB*dBB * this->L0;

            // Interpolate q and s to the face
            real_t q_face = 0, s_face = 0, wsum = 0;
            if (ir < nr) {
                const real_t w = deltaRadialFlux[ir];
                q_face += w * q_c[ir];
                s_face += w * s_c[ir];
                wsum += w;
            }
            if (ir > 0) {
                const real_t w = 1.0 - deltaRadialFlux[ir];
                q_face += w * q_c[ir-1];
                s_face += w * s_c[ir-1];
                wsum += w;
            }
            if (wsum > 0) { q_face /= wsum; s_face /= wsum; }

            // L_S = R*q/s (PoP 2024)
            const real_t s_min = 1e-6;
            real_t s_use = s_face;
            if (!std::isfinite(s_use) || std::fabs(s_use) < s_min)
                s_use = (s_use >= 0 ? s_min : -s_min);

            const real_t LS = (std::isfinite(q_face) ? (R0 * q_face / s_use) : 0.0);

            // LK = LS * (kh^2 * DFL * LS / 3)^(-1/3) (PoP 2024, QL limit)
            if (LS > 0 && std::isfinite(LS) && kh > 0 && std::isfinite(kh) && DFL > 0 && std::isfinite(DFL)) {
                const real_t X = kh*kh * DFL * LS / 3.0;
                if (X > 0 && std::isfinite(X))
                    LK = LS * std::pow(X, -1.0/3.0);
            }

            // Omega_ce = eB/me, rho_e = vTe/Omega_ce, lambda_mfp = vTe/nuc
            const real_t Omega_ce = Constants::ec * B_FSA / Constants::me;
            const real_t rho_e = (Omega_ce > 0 && std::isfinite(Omega_ce)) ? (vTe / Omega_ce) : 0.0;
            const real_t lambda_mfp = vTe / nuc;

            // LKdelta = LK * ln[ (Omega/nu) * (lcr/LK) ]
            if (LK > 0 && std::isfinite(LK) && lcr > 0 && std::isfinite(lcr) && Omega_ce > 0 && std::isfinite(Omega_ce)) {
                real_t arg = (Omega_ce / nuc) * (lcr / LK);
                arg = std::max(arg, 1.0 + 1e-12);   // keep log() well-defined
                LKdelta = LK * std::log(arg);
                if (!std::isfinite(LKdelta) || LKdelta < 0) LKdelta = 0.0;
            }

            // eta = 1 / (1 + LKdelta/lambda_mfp + lcr/rho_e)
            real_t denom = 1.0;
            if (lambda_mfp > 0 && std::isfinite(lambda_mfp))
                denom += LKdelta / lambda_mfp;
            if (rho_e > 0 && std::isfinite(rho_e))
                denom += lcr / rho_e;

            eta = 1.0 / denom;
            if (!std::isfinite(eta)) eta = 0.0;
            eta = std::max((real_t)0.0, std::min((real_t)1.0, eta));
        }

        // (NEW) store eta-related quantities for Jacobian
        this->eta_cache[ir]     = eta;
        this->LK_cache[ir]      = LK;
        this->LKdelta_cache[ir] = LKdelta;
        // ===========================================================

        // (NEW) chiQL uses eta factor: chiQL = eta * chiQL0
        real_t chiQL = 0;
        if (Theta > 0 && std::isfinite(Theta))
            chiQL = eta * qR0 * D * n * std::sqrt(Theta) * (1.0 - 5.0/8.0*Theta);

        // Lperp model (existing user extension)
        real_t Lperp = 0;
        if (dBB > 0 && std::isfinite(dBB))
            Lperp = 0.798 * Lcr / dBB;

        // ===================== soft-min connection (NEW) =====================
        const real_t p = this->pSoft;

        // chi1 = eta*chiQL0
        real_t chi1 = (chiQL > 0 ? chiQL : 0.0);

        // chi2 = chiQL*(Lperp/L0)  (strong stochastic limit scaling)
        real_t chi2 = 0.0;
        if (chiQL > 0 && Lperp > 0 && std::isfinite(Lperp) && this->L0 > 0 && std::isfinite(this->L0))
            chi2 = chiQL * (Lperp / this->L0);

        // chi3 = (3/2)*chiQL*vTe/(nuc*L0)  (collisional reduction / finite mean-free-path)
        real_t chi3 = 0.0;
        if (chiQL > 0 && vTe > 0 && std::isfinite(vTe) && nuc > 0 && std::isfinite(nuc)
            && this->L0 > 0 && std::isfinite(this->L0))
            chi3 = 1.5 * chiQL * (vTe / (nuc * this->L0));

        std::vector<real_t> chi_list = {chi1, (chi2>0?chi2:0.0), (chi3>0?chi3:0.0)};
        const real_t chiEff = softmin_pos(chi_list, p);

        Drr(ir, 0, 0) += chiEff;

        // (NEW) cache chi branches for Jacobian evaluation
        this->chiQL_cache[ir]  = (chi1>0 ? chi1 : 0.0);
        this->chi2_cache[ir]   = (chi2>0 ? chi2 : 0.0);
        this->chi3_cache[ir]   = (chi3>0 ? chi3 : 0.0);
        this->chiEff_cache[ir] = chiEff;
        // ====================================================================
    }
}


/**
 * Set jacobian of diffusion coefficients for this diffusion term.
 *
 * derivId:    ID of the quantity with respect to which the coefficient should
 *             be differentiated.
 * nMultiples: (not used).
 */
void HeatTransportRechesterRosenbluth::SetPartialDiffusionTerm(
    len_t derivId, len_t
) {
    ResetDifferentiationCoefficients();

    const len_t nr   = this->grid->GetNr();
    const real_t mc2 = Constants::mc2inEV;

    const real_t *ncold = unknowns->GetUnknownData(this->id_n_cold);
    const real_t *Tcold = unknowns->GetUnknownData(this->id_T_cold);
    const real_t *j_tot = unknowns->GetUnknownData(this->id_j_tot);

    const real_t *tauTh = nullptr;
    if (this->REFluid != nullptr)
        tauTh = this->REFluid->GetElectronCollisionTimeThermal();

    FVM::RadialGrid *rg = this->grid->GetRadialGrid();
    const real_t R0 = rg->GetR0();

    const real_t p = this->pSoft;
    const real_t EPSCHI = 1e-60;
    const real_t EPSNU  = 1e-60;

    for (len_t ir = 0; ir < nr+1; ir++) {

        // Face interpolation of T and n
        real_t T = 0, n = 0;
        if (ir < nr) {
            T += deltaRadialFlux[ir] * Tcold[ir];
            n += deltaRadialFlux[ir] * ncold[ir];
        }
        if (ir > 0) {
            T += (1.0 - deltaRadialFlux[ir]) * Tcold[ir-1];
            n += (1.0 - deltaRadialFlux[ir]) * ncold[ir-1];
        }

        const real_t Theta = T / mc2;
        if (!(Theta > 0) || !std::isfinite(Theta)) {
            dDrr(ir,0,0) = 0;
            continue;
        }
        const real_t sqrtTheta = std::sqrt(Theta);

        // qR0 (treated constant wrt n,T in this Jacobian)
        real_t qR0 = 1.0;
        if (j_tot != nullptr) {
            len_t ir_q = ir;
            if (ir_q >= nr && nr > 0) ir_q = nr-1;

            const real_t mu0Ip =
                Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir_q, rg, j_tot);
            qR0 = rg->SafetyFactorNormalized(ir_q, mu0Ip);
        }

        // nuc at face
        real_t nuc = 0;
        if (tauTh != nullptr) {
            real_t tau_face = 0, wsum = 0;

            if (ir < nr && std::isfinite(tauTh[ir])) {
                const real_t w = deltaRadialFlux[ir];
                tau_face += w * tauTh[ir];
                wsum += w;
            }
            if (ir > 0 && std::isfinite(tauTh[ir-1])) {
                const real_t w = 1.0 - deltaRadialFlux[ir];
                tau_face += w * tauTh[ir-1];
                wsum += w;
            }

            if (wsum > 0 && tau_face > 0 && std::isfinite(tau_face))
                nuc = 1.0 / tau_face;
        }

        const real_t vTe = Constants::c * std::sqrt(2.0*Theta);

        // Cached chi branches and chiEff from Rebuild()
        const real_t chi1   = this->chiQL_cache[ir];
        const real_t chi2   = this->chi2_cache[ir];
        const real_t chi3   = this->chi3_cache[ir];
        const real_t chiEff = this->chiEff_cache[ir];

        if (!(chiEff > 0) || !std::isfinite(chiEff)) {
            dDrr(ir,0,0) = 0;
            continue;
        }

        // Cached eta-model quantities (NEW)
        const real_t eta     = this->eta_cache[ir];
        const real_t LK      = this->LK_cache[ir];
        const real_t LKdelta = this->LKdelta_cache[ir];

        // chiQL0 = qR0 * D * n * sqrt(Theta) * (1 - 5/8 Theta)
        // NOTE: D is stored in this->dD[ir].
        real_t chiQL0 = qR0 * this->dD[ir] * n * sqrtTheta * (1.0 - 5.0/8.0*Theta);
        if (!(chiQL0 > 0) || !std::isfinite(chiQL0)) chiQL0 = 0;

        // Derivatives of collisionless chiQL0
        const real_t dchiQL0_dn =
            qR0 * this->dD[ir] * sqrtTheta * (1.0 - 5.0/8.0*Theta);

        const real_t dchiQL0_dT =
            qR0 * this->dD[ir]/mc2 * n
            * (0.5 - 15.0/16.0*Theta) / sqrtTheta;

        // ===================== deta/dX (NEW, partial) =====================
        // We include only eta dependence via nuc and vTe:
        //   Den = 1 + LKdelta*nuc/vTe + lcr*Omega_ce/vTe
        //   eta = 1/Den
        // and we ignore derivatives via current (q,s -> LK) and geometry.
        real_t deta_dX = 0.0;
        if (eta > 0 && std::isfinite(eta) && nuc > 0 && std::isfinite(nuc) && vTe > 0 && std::isfinite(vTe)
            && T > 0 && std::isfinite(T) && n > 0 && std::isfinite(n)) {

            // Keep lcr and Omega_ce fixed wrt n,T in this Jacobian.
            const real_t a_minor = rg->GetR_f(nr);
            const real_t m_pol   = 30.0;
            const real_t lcr = (m_pol > 0 && std::isfinite(m_pol)) ? (a_minor / m_pol) : 0.0;

            const real_t B_FSA = rg->GetFSA_B_f(ir);
            const real_t Omega_ce = Constants::ec * B_FSA / Constants::me;

            const real_t term1 = (LKdelta > 0 && std::isfinite(LKdelta)) ? (LKdelta * nuc / vTe) : 0.0;
            const real_t term2 = (lcr > 0 && std::isfinite(lcr) && Omega_ce > 0 && std::isfinite(Omega_ce))
                ? (lcr * Omega_ce / vTe) : 0.0;

            // Spitzer-like scaling used consistently with existing code:
            //   nuc ∝ n T^{-3/2}  => dnuc/dn = nuc/n,  dnuc/dT = -1.5*nuc/T
            const real_t dnuc_dn = nuc / n;
            const real_t dnuc_dT = -1.5 * nuc / T;

            // vTe = const * sqrt(T) => dvTe/dT = vTe/(2T)
            const real_t dvTe_dT = 0.5 * vTe / T;

            // dDen/dn = ((LKdelta - LK)/vTe) * (dnuc/dn)
            const real_t dDen_dn = ((LKdelta - LK) / vTe) * dnuc_dn;

            // dDen/dT = ((LKdelta - LK)/vTe) * dnuc/dT  + (dDen/dvTe)*dvTe/dT
            // where dDen/dvTe = -(term1+term2)/vTe
            const real_t dDen_dvTe = -(term1 + term2) / vTe;
            const real_t dDen_dT = ((LKdelta - LK) / vTe) * dnuc_dT + dDen_dvTe * dvTe_dT;

            if (derivId == this->id_n_cold)
                deta_dX = -eta*eta * dDen_dn;
            else if (derivId == this->id_T_cold)
                deta_dX = -eta*eta * dDen_dT;
        }
        // ================================================================

        // dchiQL/dX with chiQL = eta*chiQL0:
        //   dchiQL/dX = eta*dchiQL0/dX + chiQL0*deta/dX
        real_t dchiQL_dX = 0.0;
        if (derivId == this->id_n_cold)
            dchiQL_dX = eta * dchiQL0_dn + chiQL0 * deta_dX;
        else if (derivId == this->id_T_cold)
            dchiQL_dX = eta * dchiQL0_dT + chiQL0 * deta_dX;
        else {
            dDrr(ir,0,0) = 0;
            continue;
        }

        // soft-min partials wrt each chi_i
        const real_t dEff_dChi1 = dsoftmin_dxi(chiEff, chi1, p);
        const real_t dEff_dChi2 = dsoftmin_dxi(chiEff, chi2, p);
        const real_t dEff_dChi3 = dsoftmin_dxi(chiEff, chi3, p);

        // Build dchi1, dchi2, dchi3 (NEW: consistent with eta-corrected chiQL)
        // chi1 = chiQL
        const real_t dchi1 = dchiQL_dX;

        // chi2 = chiQL*(Lperp/L0), with d(Lperp)/d(n,T)=0 in this model
        real_t dchi2 = 0.0;
        if (chi2 > 0 && std::isfinite(chi2))
            dchi2 = (chi2 / std::max(chi1, EPSCHI)) * dchiQL_dX;

        // chi3 = (3/2)*chiQL*vTe/(nuc*L0)
        real_t dchi3 = 0.0;
        if (chi3 > 0 && std::isfinite(chi3) && nuc > 0 && std::isfinite(nuc) && vTe > 0 && std::isfinite(vTe)
            && this->L0 > 0 && std::isfinite(this->L0)) {

            // dvTe/dX: only T-derivative
            real_t dvTe_dX = 0.0;
            if (derivId == this->id_T_cold && T > 0 && std::isfinite(T))
                dvTe_dX = vTe / (2.0*T);

            // dnuc/dX (consistent with Spitzer scaling used above)
            real_t dnuc_dX = 0.0;
            if (derivId == this->id_n_cold) {
                if (n > 0 && std::isfinite(n))
                    dnuc_dX = nuc / n;
            } else if (derivId == this->id_T_cold) {
                if (T > 0 && std::isfinite(T))
                    dnuc_dX = -1.5 * nuc / T;
            }

            const real_t invNuL0 = 1.0 / (std::max(nuc, EPSNU) * this->L0);
            const real_t C = 1.5;

            // chi3 = C * chiQL * vTe * invNuL0
            // dchi3 = C*(vTe*invNuL0)*dchiQL + C*(chiQL*invNuL0)*dvTe - C*(chiQL*vTe*invNuL0/nuc)*dnuc
            dchi3 =
                C * (vTe * invNuL0) * dchiQL_dX
              + C * (std::max(chi1, (real_t)0.0) * invNuL0) * dvTe_dX
              - C * (std::max(chi1, (real_t)0.0) * vTe * invNuL0 / std::max(nuc, EPSNU)) * dnuc_dX;
        }

        // Final chain rule
        dDrr(ir,0,0) = dEff_dChi1*dchi1 + dEff_dChi2*dchi2 + dEff_dChi3*dchi3;
    }
}
