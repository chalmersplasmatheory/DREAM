/**
 * Implementation of equation term to the heating of cold electrons and ions from an NBI beam
 */
#include <cmath>
#include "DREAM/NBIHandler.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "DREAM/ADAS.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Interpolator1D.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Constants.hpp"
#include <tuple>
#include <vector>


using namespace DREAM;
/**
 * Constructor initializing member variables
 */
NBIHandler::NBIHandler(FVM::Grid* grid, ADAS* adas, IonHandler* ions
): grid(grid), radialGrid(grid->GetRadialGrid()), adas(adas), ions(ions) {

    nr = radialGrid->GetNr();
    nZ = ions->GetNZ();
    len_t Zmax_max = 0;
    for (len_t iz = 0; iz < nZ; ++iz)
        Zmax_max = std::max(Zmax_max, ions->GetZ(iz));
    nCharge = Zmax_max + 1;
    size_t maxidx =
        static_cast<size_t>(nr) * static_cast<size_t>(nZ) * static_cast<size_t>(nCharge);

    NBIHeatTerm_e = new real_t[nr];
    NBIHeatTerm_i = new real_t[nr];
    Deposition_profile = new real_t[nr];
    Deposition_profile_times_Vprime = new real_t[nr];
    dV_beam_prime_tot = new real_t[nr];
    H_r_dTe = new real_t[nr];
    H_r_dne = new real_t[nr];
    d_NBIHeatTerm_e_d_Te = new real_t[nr];
    d_NBIHeatTerm_e_d_ne = new real_t[nr];
    d_NBIHeatTerm_i_d_Te = new real_t[nr];
    d_NBIHeatTerm_i_d_ne = new real_t[nr];
    d_NBIHeatTerm_i_d_T_ij = new real_t[maxidx];
    d_NBIHeatTerm_i_d_n_ij = new real_t[maxidx];
    d_NBIHeatTerm_e_d_T_ij = new real_t[maxidx];
    d_NBIHeatTerm_e_d_n_ij = new real_t[maxidx];
    H_r_dn_ij = new real_t[maxidx];
    H_r_dT_ij = new real_t[maxidx];
}

/**
 * Configure the NBIHandler from NBI settings. Called at the start of the simulation.
 */
void NBIHandler::ConfigureFromSettings(
    Settings* , FVM::UnknownQuantityHandler* unknowns, real_t s_max, real_t r_beam,
    const real_t P0[3], const real_t n[3], const real_t energy_fractions[3],
    real_t Ti_beam, real_t m_i_beam, FVM::Interpolator1D* j_B_profile,
    real_t R0,  int gaussian_profile, FVM::Interpolator1D* Power_Profile
) {

    this->unknowns = unknowns;
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);

    // Set integration step size
    this->ntheta = 100;
    this->dtheta = 2.0 * M_PI / ntheta;

    // Set beam parameters
    this->s_max = s_max;
    this->r_beam = r_beam;
    this->Ti_beam = Ti_beam;
    this->m_i_beam = m_i_beam;
    this->R0 = R0;
    if (P0 == nullptr)
        throw DREAMException("Starting point of the beam is not set correctly");
    if (n == nullptr)
        throw DREAMException("Direction vector of the beam is not set correctly");

    for (int i = 0; i < 3; ++i) {
        this->P0[i] = P0[i];
        this->n[i] = n[i];
        this->energy_fractions[i] = energy_fractions[i];
    }

    // Compute the beam basis vectors e1 and e2 and the start and stop of the beam in beam coordinates, initialized to extreme values
    this->s_start = 1e10;
    this->s_stop = 0;
    PrecomputeBeamBasisVectors();
    PrecomputeFluxSurfaces();

    if (s_start > s_stop) 
        throw DREAMException("Beam does not intersect with any flux surface.");
    if (s_max < s_stop) 
        throw DREAMException("Beam does not extend all the way through the plasma.");
    if (s_start < 0) 
        throw DREAMException("Beam starts in the plasma (s_start < 0).");

    this->j_B_profile = j_B_profile;
    this->Power_Profile = Power_Profile;
    this->gaussian_profile = gaussian_profile;

    if (Power_Profile == nullptr) 
        throw DREAMException("NBI: Power profile must be provided (P_NBI).");

    if (Power_Profile->GetNx() == 0) 
        throw DREAMException("NBI: Power profile is empty.");

    this->n_beam_radius = 25;
    this->d_beam_radius = r_beam / n_beam_radius;
    this->n_beam_theta = 25;
    this->d_beam_theta = 2.0 * M_PI / n_beam_theta;
    this->n_beam_s = 50;
    this->d_beam_s = (s_stop - s_start) / n_beam_s;

    //Precompute the ir-flux surface mapping
    PrecomputeBeamMapLUT();

    // Trapezoidal integration to get total beam current, if a jB_profile is set
    this->I_B = 0.0;
    const real_t* r_beam_grid = this->j_B_profile->GetX();
    len_t N = this->j_B_profile->GetNx();
    for (len_t i = 0; i < N - 1; ++i) {
        real_t r1 = r_beam_grid[i];
        real_t r2 = r_beam_grid[i + 1];
        real_t dr = r2 - r1;

        // Evaluate current density at physical r values
        real_t j1 = this->j_B_profile->Eval(r1)[0];
        real_t j2 = this->j_B_profile->Eval(r2)[0];

        // Trapezoidal integration
        this->I_B += 2.0 * M_PI * 0.5 * (j1 * r1 + j2 * r2) * dr;
    }
}

/**
 * Calculate the beam current density j_B divided by I_B.
 */
real_t NBIHandler::Calculate_jB_IB(real_t r_B, real_t theta_B) {
    if (this->gaussian_profile == OptionConstants::EQTERM_NBI_GAUSSIAN_PROFILE_ITER) {
        // ITER-like beam cross-section, elliptical Gaussian profile
        real_t Dx = 1.0;
        real_t Dy = 0.5;
        real_t A = 4.0 / (M_PI * Dx * Dy);
        real_t cs = std::cos(theta_B);
        real_t sn = std::sin(theta_B);
        real_t exponent = -4.0 * r_B * r_B * (cs * cs / (Dx * Dx) + sn * sn / (Dy * Dy));
        return A * std::exp(exponent);

    } else if (this->gaussian_profile == OptionConstants::EQTERM_NBI_GAUSSIAN_PROFILE_TCV) {
        // TCV-like beam cross-section, elliptical Gaussian profile
        real_t Dx = 0.3;
        real_t Dy = 0.1;
        real_t A = 4 / (M_PI * Dx * Dy);
        real_t cs = std::cos(theta_B);
        real_t sn = std::sin(theta_B);
        real_t exponent = -4.0 * r_B * r_B * (cs * cs / (Dx * Dx) + sn * sn / (Dy * Dy));
        real_t jB_divided_IB = A * std::exp(exponent);
        return jB_divided_IB;
    } else if (this->gaussian_profile == OptionConstants::EQTERM_NBI_GAUSSIAN_PROFILE_CUSTOM) {
        // Default case: use jB_profile
        const real_t* j_B_values = j_B_profile->Eval(r_B);
        real_t j_B_r = j_B_values[0];
        return j_B_r / this->I_B;
    }
    else {
        throw DREAMException("NBI: Unrecognized Gaussian profile type %d.", this->gaussian_profile);
    }
}

/**
 * Precomputes the flux surfaces in R and Z. Also find the start and end of the beam in beam coordinates.
 */
void NBIHandler::PrecomputeFluxSurfaces() {

    cachedFluxSurfaces.resize(nr);
    for (len_t ir = 0; ir < nr; ++ir) {
        cachedFluxSurfaces[ir].resize(ntheta);
        for (real_t itheta = 0; itheta < ntheta; ++itheta) {
            real_t theta = itheta * dtheta;
            real_t R = R0 + radialGrid->GetFluxSurfaceRMinusR0_theta(ir, theta);
            real_t Z = radialGrid->GetFluxSurfaceZMinusZ0_theta(ir, theta);
            cachedFluxSurfaces[ir][itheta] = {R, Z};

            real_t nphi = 200;
            real_t dphi = 2.0 * M_PI / nphi;
            if (ir == nr - 1) {
                for (real_t iphi = 0; iphi < nphi; ++iphi) {
                    real_t phi = iphi * dphi;
                    real_t x = R * std::cos(phi);
                    real_t y = R * std::sin(phi);
                    real_t z = Z;
                    real_t r_B, theta_B, s_B;
                    CartesianToCylindrical(x, y, z, r_B, theta_B, s_B);
                    if (r_B < r_beam) {
                        if (s_B > s_stop)
                            s_stop = s_B;
                        if (s_B < s_start)
                            s_start = s_B;
                    }
                }
            }
        }
    }
}

/**
 * Precomputes the beam basis vectors of beam e1 and e2.
 */
void NBIHandler::PrecomputeBeamBasisVectors() {
    real_t norm_n = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (int i = 0; i < 3; ++i)
        n_norm[i] = n[i] / norm_n;

    std::array<real_t, 3> a;
    if (std::abs(n_norm[2]) < 0.9)
        a = {0, 0, 1};
    else
        a = {1, 0, 0};
    real_t a_dot_n = a[0] * n_norm[0] + a[1] * n_norm[1] + a[2] * n_norm[2];
    for (int i = 0; i < 3; ++i)
        e1[i] = a[i] - a_dot_n * n_norm[i];

    real_t norm_e1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (int i = 0; i < 3; ++i)
        e1[i] /= norm_e1;

    e2[0] = n_norm[1] * e1[2] - n_norm[2] * e1[1];
    e2[1] = n_norm[2] * e1[0] - n_norm[0] * e1[2];
    e2[2] = n_norm[0] * e1[1] - n_norm[1] * e1[0];
}

/**
 * Precomputes discretized points in the beam and maps the intersections of the flux surfaces.
 */
void NBIHandler::PrecomputeBeamMapLUT() {
    ir_lut.resize(static_cast<size_t>(n_beam_theta) * n_beam_radius * n_beam_s,
                  static_cast<len_t>(-1));

    for (len_t it = 0; it < n_beam_theta; ++it) {
        for (len_t irad = 0; irad < n_beam_radius; ++irad) {
            for (len_t is = 0; is < n_beam_s; ++is) {
                real_t beam_theta = it * d_beam_theta;
                real_t beam_radius = (irad) *d_beam_radius;
                real_t beam_s = s_start + is * d_beam_s;

                len_t ir_now =
                    CalculatePencilBeamFindFlux(beam_s, beam_radius, beam_theta);
                ir_lut[lut_index(it, irad, is)] = ir_now;
            }
        }
    }
}

/**
 * Destructor
 */
NBIHandler::~NBIHandler() {
    delete[] NBIHeatTerm_e;
    delete[] NBIHeatTerm_i;
    delete[] Deposition_profile;
    delete[] Deposition_profile_times_Vprime;
    delete[] dV_beam_prime_tot;
    delete[] H_r_dTe;
    delete[] H_r_dne;
    delete[] d_NBIHeatTerm_e_d_Te;
    delete[] d_NBIHeatTerm_e_d_ne;
    delete[] d_NBIHeatTerm_i_d_Te;
    delete[] d_NBIHeatTerm_i_d_ne;
    delete[] d_NBIHeatTerm_i_d_T_ij;
    delete[] d_NBIHeatTerm_i_d_n_ij;
    delete[] d_NBIHeatTerm_e_d_T_ij;
    delete[] d_NBIHeatTerm_e_d_n_ij;
    delete[] H_r_dn_ij;
    delete[] H_r_dT_ij;
}

/**
 * Build: Called at the start of every timestep.
 * Compute added heating term for each timestep and flux surface radius
 */
void NBIHandler::Build(const real_t t, const real_t,
    FVM::UnknownQuantityHandler* unknowns, real_t Ti_beam_input
) {

    // Setting values to zero since build is called multiple times per timestep
    const real_t* VpVol = this->radialGrid->GetVpVol();
    const real_t* dr = this->radialGrid->GetDr();
    for (len_t ir = 0; ir < nr; ++ir) {
        NBIHeatTerm_e[ir] = 0.0;
        NBIHeatTerm_i[ir] = 0.0;
        Deposition_profile[ir] = 0.0;
        Deposition_profile_times_Vprime[ir] = 0.0;
        dV_beam_prime_tot[ir] = 0.0;
        H_r_dTe[ir] = 0.0;
        H_r_dne[ir] = 0.0;
        d_NBIHeatTerm_e_d_Te[ir] = 0.0;
        d_NBIHeatTerm_e_d_ne[ir] = 0.0;
        d_NBIHeatTerm_i_d_Te[ir] = 0.0;
        d_NBIHeatTerm_i_d_ne[ir] = 0.0;
        for (len_t iz = 0; iz < nZ; iz++) {
            for (len_t Z0 = 0; Z0 <= ions->GetZ(iz); Z0++) {
                H_r_dT_ij[idx(ir, iz, Z0)] = 0.0;
                H_r_dn_ij[idx(ir, iz, Z0)] = 0.0;
                d_NBIHeatTerm_i_d_n_ij[idx(ir, iz, Z0)] = 0.0;
                d_NBIHeatTerm_i_d_T_ij[idx(ir, iz, Z0)] = 0.0;
                d_NBIHeatTerm_e_d_n_ij[idx(ir, iz, Z0)] = 0.0;
                d_NBIHeatTerm_e_d_T_ij[idx(ir, iz, Z0)] = 0.0;
            }
        }
    }
    const real_t* p = this->Power_Profile->Eval(t);
    real_t powerNow = p[0];
    
    // If the power is not 0, then compute the calculation
    if (powerNow != 0) {

        // Compute the deposition profile
        ComputeDepositionProfile(unknowns, Ti_beam_input);

        // Compute the deposited fraction of the beam
        real_t frac = 0.0;
        for (len_t ir = 1; ir < nr; ir++) {
            frac += Deposition_profile_times_Vprime[ir] * d_beam_radius;
        }
        real_t fi, fe, dfe_dne, dfe_dTe, dfi_dne, dfi_dTe;
        std::vector<real_t> dfe_dn_ij(nZ), dfi_dn_ij(nZ), dfe_dT_ij(nZ), dfi_dT_ij(nZ);
        for (len_t ir = 0; ir < nr; ++ir) {

            // Calculate ion and electron fractions and their derivatives
            IonElectronFractions(unknowns, ir, fi, fe, dfe_dne, dfe_dTe, dfe_dn_ij,
                                 dfe_dT_ij, dfi_dne, dfi_dTe, dfi_dn_ij, dfi_dT_ij,
                                 Ti_beam_input);
            real_t K = powerNow * d_beam_radius / (VpVol[ir] * dr[ir]);

            // Set heating terms
            NBIHeatTerm_e[ir] = fe * powerNow * Deposition_profile_times_Vprime[ir] * d_beam_radius / (VpVol[ir] * dr[ir]);
            NBIHeatTerm_i[ir] = fi * powerNow * Deposition_profile_times_Vprime[ir] * d_beam_radius / (VpVol[ir] * dr[ir]);

            //Set derivatives of heating terms
            d_NBIHeatTerm_e_d_Te[ir] = K * (fe * H_r_dTe[ir] + dfe_dTe * Deposition_profile_times_Vprime[ir]);
            d_NBIHeatTerm_e_d_ne[ir] = K * (fe * H_r_dne[ir] + dfe_dne * Deposition_profile_times_Vprime[ir]);
            d_NBIHeatTerm_i_d_Te[ir] = K * (fi * H_r_dTe[ir] + dfi_dTe * Deposition_profile_times_Vprime[ir]);
            d_NBIHeatTerm_i_d_ne[ir] = K * (fi * H_r_dne[ir] + dfi_dne * Deposition_profile_times_Vprime[ir]);
            for (len_t iz = 0; iz < nZ; iz++) {
                for (len_t Z0 = 0; Z0 <= ions->GetZ(iz); Z0++) {
                    d_NBIHeatTerm_i_d_n_ij[idx(ir, iz, Z0)] = K * (fi * H_r_dn_ij[idx(ir, iz, Z0)] + dfi_dn_ij[iz] * Deposition_profile_times_Vprime[ir]);
                    d_NBIHeatTerm_i_d_T_ij[idx(ir, iz, Z0)] = K * (fi * H_r_dT_ij[idx(ir, iz, Z0)] + dfi_dT_ij[iz] * Deposition_profile_times_Vprime[ir]);
                    d_NBIHeatTerm_e_d_n_ij[idx(ir, iz, Z0)] = K * (fe * H_r_dn_ij[idx(ir, iz, Z0)] + dfe_dn_ij[iz] * Deposition_profile_times_Vprime[ir]);
                    d_NBIHeatTerm_e_d_T_ij[idx(ir, iz, Z0)] = K * (fe * H_r_dT_ij[idx(ir, iz, Z0)] + dfe_dT_ij[iz] * Deposition_profile_times_Vprime[ir]);
                }
            }
        }

    } else {
        // If power is 0, set everything to zero
        for (len_t ir = 0; ir < nr; ir++) {
            NBIHeatTerm_e[ir] = 0;
            NBIHeatTerm_i[ir] = 0;
            H_r_dTe[ir] = 0;
            H_r_dne[ir] = 0;
            d_NBIHeatTerm_e_d_Te[ir] = 0;
            d_NBIHeatTerm_e_d_ne[ir] = 0;
            d_NBIHeatTerm_i_d_Te[ir] = 0;
            d_NBIHeatTerm_i_d_ne[ir] = 0;
            for (len_t iz = 0; iz < nZ; iz++) {
                for (len_t Z0 = 0; Z0 <= ions->GetZ(iz); Z0++) {
                    H_r_dT_ij[idx(ir, iz, Z0)] = 0;
                    H_r_dn_ij[idx(ir, iz, Z0)] = 0;
                    d_NBIHeatTerm_i_d_n_ij[idx(ir, iz, Z0)] = 0;
                    d_NBIHeatTerm_i_d_T_ij[idx(ir, iz, Z0)] = 0;
                    d_NBIHeatTerm_e_d_n_ij[idx(ir, iz, Z0)] = 0;
                    d_NBIHeatTerm_e_d_T_ij[idx(ir, iz, Z0)] = 0;
                }
            }
        }
    }
}

/**
 * Computes and mean free path variable for the NBI heating term
 */
void NBIHandler::ComputeMeanFreePath(
    len_t ir, real_t ncold, real_t Tcold,
    real_t& lambda_s, real_t& dlambda_dI_e,
    std::vector<real_t>& dlambda_dI_ij,
    real_t& dlambda_dne, std::vector<real_t>& dlambda_dn_ij,
    std::vector<real_t>& dI_ij_dT_ij, real_t& dI_e_dTe,
    real_t Ti_beam_input
) {
    real_t Z0 = 0; // Neutral beam, so Z0 = 0
    real_t Zion = 1; // Deterium beam
    ADASRateInterpolator* scd = adas->GetSCD(Zion);
    dI_e_dTe = scd->Eval_deriv_T(Z0, ncold, Tcold);

    real_t dI_e_dne = scd->Eval_deriv_n(Z0, ncold, Tcold);
    real_t I_e = scd->Eval(Z0, ncold, Tcold);
    real_t v_NBI = std::sqrt(2.0 * Ti_beam_input / m_i_beam);
    real_t total_CX = 0.0;

    // Loop over all atomic species
    for (len_t iz = 0; iz < nZ; iz++) {
        len_t Zmax = ions->GetZ(iz);

        real_t niZ_tot = 0.0;
        for (len_t Z0 = 0; Z0 <= Zmax; Z0++) {
            niZ_tot += ions->GetIonDensity(ir, iz, Z0);
        }
        // Calculate temperature for individual species
        real_t TiZ = (unknowns->GetUnknownData(id_ion_temperature)[iz * nr + ir]) /(1.5 * niZ_tot * Constants::ec);

        // Calculate total CX rate
        ADASRateInterpolator* ccd = adas->GetCCD(Zmax);
        for (len_t Z0 = 0; Z0 <= Zmax; Z0++) {
            real_t niZ = ions->GetIonDensity(ir, iz, Z0);

            // Evaluate CX rate coefficient for every species and charge state
            real_t I_ij = ccd->Eval(Z0, niZ, TiZ);
            total_CX += I_ij * niZ;
        }
    }

    // Do the same calculation for the derivatives where the total_CX is needed
    for (len_t iz = 0; iz < nZ; iz++) {
        len_t Zmax = ions->GetZ(iz);

        ADASRateInterpolator* ccd = adas->GetCCD(Zmax);
        real_t niZ_tot = 0.0;
        for (len_t Z0 = 0; Z0 <= Zmax; Z0++) {
            niZ_tot += ions->GetIonDensity(ir, iz, Z0);
        }

        real_t TiZ = (unknowns->GetUnknownData(id_ion_temperature)[iz * nr + ir]) /
            (1.5 * niZ_tot * Constants::ec);
        for (len_t Z0 = 0; Z0 <= Zmax; Z0++) {
            real_t niZ = ions->GetIonDensity(ir, iz, Z0);
            real_t I_ij = ccd->Eval(Z0, niZ, TiZ);
            real_t dI_ij_dn_ij = ccd->Eval_deriv_n(Z0, niZ, TiZ);
            
            dI_ij_dT_ij[idx(ir, iz, Z0)] = ccd->Eval_deriv_T(Z0, niZ, TiZ);
            dlambda_dI_ij[idx(ir, iz, Z0)] = -v_NBI * niZ / ((ncold * I_e + total_CX) * (ncold * I_e + total_CX));
            dlambda_dn_ij[idx(ir, iz, Z0)] = -v_NBI * (I_ij + dI_ij_dn_ij * niZ) / ((ncold * I_e + total_CX) * (ncold * I_e + total_CX));
        }
    }

    lambda_s = v_NBI / (ncold * (I_e) + total_CX);
    dlambda_dI_e = -v_NBI * ncold / ((ncold * I_e + total_CX) * (ncold * I_e + total_CX));
    dlambda_dne = -v_NBI * (I_e + dI_e_dne * ncold) / ((ncold * I_e + total_CX) * (ncold * I_e + total_CX));
}

/**
 * Computes the flux surface index for a given pencil beam position (s_B, r_B, theta_B).
 * This function checks if the pencil beam intersects with any flux surface and returns the index of the first one it finds.
 */
int_t NBIHandler::CalculatePencilBeamFindFlux(real_t s_B, real_t r_B, real_t theta_B) {
    // Convert the beam coordinates (s_B, r_B, theta_B) to Cartesian coordinates
    real_t dx = r_B * cos(theta_B);
    real_t dy = r_B * sin(theta_B);
    std::array<real_t, 3> P0_current = {P0[0] + dx * e1[0] + dy * e2[0],
                                        P0[1] + dx * e1[1] + dy * e2[1],
                                        P0[2] + dx * e1[2] + dy * e2[2]};

    int_t best_ir = -1;
    real_t smallest_diff = std::numeric_limits<real_t>::infinity();

    for (len_t ir = 1; ir < nr; ir++) {

        for (len_t i = 0; i < ntheta; i++) {
            len_t j = (i + 1) % ntheta;

            const FluxSurfacePoint& pt1 = cachedFluxSurfaces[ir][i];
            const FluxSurfacePoint& pt2 = cachedFluxSurfaces[ir][j];

            real_t R1 = pt1.R, Z1 = pt1.Z;
            real_t R2 = pt2.R, Z2 = pt2.Z;

            real_t dZ = Z2 - Z1;
            real_t dR = R2 - R1;

            real_t ratio = dR / dZ;

            real_t a0 = R1 * R1 - P0_current[0] * P0_current[0] -
                P0_current[1] * P0_current[1] + 2 * R1 * (P0_current[2] - Z1) * ratio +
                (P0_current[2] - Z1) * (P0_current[2] - Z1) * ratio * ratio;

            real_t a1 = 2 * R1 * n_norm[2] * ratio +
                2 * n_norm[2] * (P0_current[2] - Z1) * ratio * ratio -
                2 * (P0_current[0] * n_norm[0] + P0_current[1] * n_norm[1]);

            real_t a2 = n_norm[2] * n_norm[2] * ratio * ratio - n_norm[0] * n_norm[0] -
                n_norm[1] * n_norm[1];

            real_t discriminant = a1 * a1 / (4 * a2 * a2) - a0 / a2;
            if (discriminant <= 0)
                continue;


            real_t sqrt_disc = std::sqrt(discriminant);
            real_t l1 = (-a1 / (2 * a2) + sqrt_disc);
            real_t l2 = (-a1 / (2 * a2) - sqrt_disc);

            auto check_solution = [&](real_t l) -> std::pair<bool, real_t> {
                if (l < 0)
                    return {false, 0};
                real_t t = (P0_current[2] - Z1 + l * n_norm[2]) / dZ;
                if (t < 0 || t > 1)
                    return {false, 0};
                return {true, std::abs(l - s_B)};
            };

            // Check both solutions and update if better match found
            auto [valid1, diff1] = check_solution(l1);

            if (valid1 && diff1 < smallest_diff) {
                smallest_diff = diff1;
                best_ir = ir;
            }

            auto [valid2, diff2] = check_solution(l2);
            if (valid2 && diff2 < smallest_diff) {
                smallest_diff = diff2;
                best_ir = ir;
            }
        }
    }
    return best_ir;
}
/**
 * Takes in the beam paramaters and turns a point on a flux surface (in x y z) into the beam coordinate system
 */
void NBIHandler::CartesianToCylindrical(real_t x, real_t y, real_t z, real_t& r,
                                        real_t& theta, real_t& s) {

    std::array<real_t, 3> d = {x - P0[0], y - P0[1], z - P0[2]};
    s = d[0] * n_norm[0] + d[1] * n_norm[1] + d[2] * n_norm[2];

    std::array<real_t, 3> d_trans;
    for (int i = 0; i < 3; ++i)
        d_trans[i] = d[i] - s * n_norm[i];

    real_t r_cos_theta = d_trans[0] * e1[0] + d_trans[1] * e1[1] + d_trans[2] * e1[2];
    real_t r_sin_theta = d_trans[0] * e2[0] + d_trans[1] * e2[1] + d_trans[2] * e2[2];
    r = std::sqrt(r_cos_theta * r_cos_theta + r_sin_theta * r_sin_theta);
    theta = atan2(r_sin_theta, r_cos_theta);
    if (theta < 0)
        theta += 2 * M_PI;
}

/**
 * Computes the deposition profile for a flux surface
 */
void NBIHandler::ComputeDepositionProfile(FVM::UnknownQuantityHandler* unknowns, real_t Ti_beam_input) {

    for (len_t ir = 0; ir < nr; ++ir) {
        Deposition_profile_times_Vprime[ir] = 0.0;
        Deposition_profile[ir] = 0.0;
        dV_beam_prime_tot[ir] = 0.0;
        H_r_dTe[ir] = 0.0;
        H_r_dne[ir] = 0.0;
        for (len_t iz = 0; iz < nZ; iz++) {
            for (len_t Z0 = 0; Z0 <= ions->GetZ(iz); Z0++) {
                H_r_dT_ij[idx(ir, iz, Z0)] = 0.0;
                H_r_dn_ij[idx(ir, iz, Z0)] = 0.0;
            }
        }
    }
    std::vector<real_t> lambda_cache(nr, -1.0); // -1 = not yet computed
    std::vector<real_t> dlambda_dI_cache(nr, -1.0);
    std::vector<real_t> dlambda_dne_cache(nr, -1.0);
    std::vector<real_t> dI_e_dTe_cache(nr, -1.0);
    std::vector<real_t> dlambda_dI_ij_cache(nr * nZ * nCharge, -1.0);
    std::vector<real_t> dlambda_dn_ij_cache(nr * nZ * nCharge, -1.0);
    std::vector<real_t> dI_ij_dT_ij_cache(nr * nZ * nCharge, -1.0);


    real_t lambda_s, dlambda_dI_e, dlambda_dne, dI_e_dTe;
    std::vector<real_t> dlambda_dI_ij(nr * nZ * nCharge, 0.0);
    std::vector<real_t> dlambda_dn_ij(nr * nZ * nCharge, 0.0);
    std::vector<real_t> dI_ij_dT_ij(nr * nZ * nCharge, 0.0);

    for (real_t i_beam_theta = 0; i_beam_theta < n_beam_theta; i_beam_theta++) {
        real_t beam_theta = i_beam_theta * d_beam_theta;
        for (real_t i_beam_radius = 0; i_beam_radius < n_beam_radius; i_beam_radius++) {

            real_t beam_radius = (i_beam_radius) *d_beam_radius;
            real_t I_s = 0.0;
            real_t I_s_squared = 0.0;
            real_t jB_divided_IB = Calculate_jB_IB(beam_radius, beam_theta);
            for (real_t i_beam_s = 0; i_beam_s < n_beam_s; i_beam_s++) {
                real_t beam_s = s_start + i_beam_s * d_beam_s;
                len_t ir_now = ir_lut[lut_index(i_beam_theta, i_beam_radius, i_beam_s)];

                if (ir_now >= nr) {
                    throw DREAMException(
                        "NBIElectronHeatTerm: ir_now out of bounds: ir_now = "
                        LEN_T_PRINTF_FMT " (s_B = %.4f, r_B = %.4f, theta_B = %.4f)",
                        ir_now, beam_s, beam_radius, beam_theta
                    );
                }
                real_t ncold = unknowns->GetUnknownData(id_ncold)[ir_now];
                real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir_now];

                // Use cached value if available
                if (lambda_cache[ir_now] < 0) {
                    ComputeMeanFreePath(
                        ir_now, ncold, Tcold, lambda_s, dlambda_dI_e,
                        dlambda_dI_ij, dlambda_dne, dlambda_dn_ij,
                        dI_ij_dT_ij, dI_e_dTe, Ti_beam_input
                    );

                    lambda_cache[ir_now] = lambda_s;
                    dlambda_dI_cache[ir_now] = dlambda_dI_e;
                    dlambda_dne_cache[ir_now] = dlambda_dne;
                    dI_e_dTe_cache[ir_now] = dI_e_dTe;

                    for (len_t iz = 0; iz < nZ; iz++) {
                        for (len_t Z0 = 0; Z0 <= ions->GetZ(iz); Z0++) {
                            dlambda_dI_ij_cache[idx(ir_now, iz, Z0)] = dlambda_dI_ij[idx(ir_now, iz, Z0)];
                            dlambda_dn_ij_cache[idx(ir_now, iz, Z0)] = dlambda_dn_ij[idx(ir_now, iz, Z0)];
                            dI_ij_dT_ij_cache[idx(ir_now, iz, Z0)] = dI_ij_dT_ij[idx(ir_now, iz, Z0)];
                        }
                    }
                } else {
                    lambda_s = lambda_cache[ir_now];
                    dlambda_dI_e = dlambda_dI_cache[ir_now];
                    dlambda_dne = dlambda_dne_cache[ir_now];
                    dI_e_dTe = dI_e_dTe_cache[ir_now];
                    for (len_t iz = 0; iz < nZ; iz++) {
                        for (len_t Z0 = 0; Z0 <= ions->GetZ(iz); Z0++) {
                            dlambda_dI_ij[idx(ir_now, iz, Z0)] = dlambda_dI_ij_cache[idx(ir_now, iz, Z0)];
                            dlambda_dn_ij[idx(ir_now, iz, Z0)] = dlambda_dn_ij_cache[idx(ir_now, iz, Z0)];
                            dI_ij_dT_ij[idx(ir_now, iz, Z0)] = dI_ij_dT_ij_cache[idx(ir_now, iz, Z0)];
                        }
                    }
                }
                I_s += 1 / lambda_s * d_beam_s;
                I_s_squared += 1 / (lambda_s * lambda_s) * d_beam_s;
                real_t survivalProb = std::exp(-I_s);

                real_t dV_beam_prime = beam_radius * d_beam_theta * d_beam_s;
                real_t H_r_theta_phi = jB_divided_IB * (1.0 / lambda_s) * survivalProb;
                Deposition_profile[ir_now] += H_r_theta_phi;

                Deposition_profile_times_Vprime[ir_now] += H_r_theta_phi * dV_beam_prime;

                dV_beam_prime_tot[ir_now] += dV_beam_prime;

                // Derivatives
                H_r_dTe[ir_now] += jB_divided_IB * dlambda_dI_e * dI_e_dTe * ((-1.0 / (lambda_s * lambda_s)) * survivalProb + 1 / lambda_s * I_s_squared * survivalProb) * dV_beam_prime;
                H_r_dne[ir_now] += jB_divided_IB * dlambda_dne * ((-1.0 / (lambda_s * lambda_s)) * survivalProb + 1 / lambda_s * I_s_squared * survivalProb) * dV_beam_prime;
                for (len_t iz = 0; iz < nZ; iz++) {
                    for (len_t Z0 = 0; Z0 <= ions->GetZ(iz); Z0++) {
                        H_r_dT_ij[idx(ir_now, iz, Z0)] += jB_divided_IB * dlambda_dI_ij[idx(ir_now, iz, Z0)] * dI_ij_dT_ij[idx(ir_now, iz, Z0)] *
                            ((-1.0 / (lambda_s * lambda_s)) * survivalProb +
                             1 / lambda_s * I_s_squared * survivalProb) *
                            dV_beam_prime;
                        H_r_dn_ij[idx(ir_now, iz, Z0)] += jB_divided_IB *
                            dlambda_dn_ij[idx(ir_now, iz, Z0)] *
                            ((-1.0 / (lambda_s * lambda_s)) * survivalProb +
                             1 / lambda_s * I_s_squared * survivalProb) *
                            dV_beam_prime;
                    }
                }
            }
        }
    }
}


/**
 *Calculate the fraction of energy that should go to the electrons and different ion species
 */
void NBIHandler::IonElectronFractions(
    FVM::UnknownQuantityHandler* unknowns, len_t ir,
    real_t& f_i, real_t& f_e, real_t& dfe_dne,
    real_t& dfe_dTe, std::vector<real_t>& dfe_dn_ij,
    std::vector<real_t>& dfe_dT_ij, real_t& dfi_dne,
    real_t& dfi_dTe, std::vector<real_t>& dfi_dn_ij,
    std::vector<real_t>& dfi_dT_ij,
    real_t Ti_beam_input
) {
    real_t m_u = Constants::mu;
    real_t evJ = Constants::ec;
    real_t me = Constants::me;

    real_t ne = unknowns->GetUnknownData(id_ncold)[ir];
    real_t Te = unknowns->GetUnknownData(id_Tcold)[ir];

    real_t Me = me / m_u;
    real_t Mb = this->m_i_beam / m_u;
    real_t ve = std::sqrt(2 * Te / Me);
    real_t vb = std::sqrt(2 * (Ti_beam_input / evJ) / Mb);

    std::vector<real_t> ni_species(nZ, 0.0);
    std::vector<real_t> Zi(nZ, 0.0);
    std::vector<real_t> Mi(nZ, 0.0);
    // Loop over all atomic species
    for (len_t iz = 0; iz < nZ; iz++) {
        len_t Zmax = ions->GetZ(iz);
        Zi[iz] = Zmax;
        Mi[iz] = ions->GetIonSpeciesMass(iz) / m_u;


        // Loop over charge states of this species
        real_t sum_ni = 0.0;
        for (len_t Z0 = 0; Z0 <= Zmax; Z0++) {
            real_t niZ = ions->GetIonDensity(ir, iz, Z0);
            sum_ni += niZ;
        }
        ni_species[iz] = sum_ni;
    }
    real_t Z1_sum = 0.0;
    for (len_t iz = 0; iz < nZ; iz++) {
        if (Mi[iz] > 0.0)
            Z1_sum += ni_species[iz] * (Zi[iz] * Zi[iz]) / Mi[iz];
    }
    real_t Z1 = (Mb / ne) * Z1_sum;
    real_t factor = (3 * std::sqrt(M_PI) * Me * Z1 / (4 * Mb));
    real_t vc = std::pow(std::max((real_t) 0, factor), (real_t) (1.0 / 3.0)) * ve;
    real_t xc = vc / vb;
    real_t x = 1 / (xc);
    real_t F_NB2 = (std::atan((2 * x - 1) / std::sqrt(3)) + M_PI / 6) / std::sqrt(3) -
        (std::log((1 + x) * (1 + x) / (1 - x + x * x))) / 6.0;

    f_i = 2 * F_NB2 / (x * x);
    f_e = 1 - f_i;

    //Derivatives calculation
    real_t F_x = 1 / std::sqrt(3) * std::atan((2 * x - 1) / std::sqrt(3)) +
        M_PI / (6 * std::sqrt(3)) - std::log((1 + x) * (1 + x) / (1 - x + x * x)) / 6;
    real_t F_x_prime = 2 / (3 + (2 * x - 1) * (2 * x - 1)) - 1 / (3 * (1 + x)) -
        (2 * x - 1) / (6 * (1 - x + x * x));
    real_t G_x = 2 / (x * x * x) * (x * F_x_prime - 2 * F_x);

    dfe_dne = -G_x * x / (2 * ne);
    dfe_dTe = G_x * x / (3 * Te);
    for (len_t iz = 0; iz < nZ; iz++) {
        dfe_dn_ij[iz] = G_x * x / 3 *
            (ions->GetZ(iz) * ions->GetZ(iz) / (ions->GetIonSpeciesMass(iz) / m_u) /
             Z1_sum);
        dfe_dT_ij[iz] = 0;
        dfi_dn_ij[iz] = -dfe_dn_ij[iz];
        dfi_dT_ij[iz] = -dfe_dT_ij[iz];
    }
    dfi_dne = -dfe_dne;
    dfi_dTe = -dfe_dTe;
}