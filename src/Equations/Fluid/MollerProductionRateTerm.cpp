#include "DREAM/Equations/Fluid/MollerProductionRateTerm.hpp"

#include <algorithm>

#include "DREAM/Settings/OptionConstants.hpp"

using namespace DREAM;

MollerProductionRateTerm::MollerProductionRateTerm(
    FVM::Grid *fluidGrid, const FVM::Grid *grid_target, const FVM::Grid *grid_primary,
    FVM::UnknownQuantityHandler *u, len_t id_f_primary, const MollerEnergyKernel *energy_kernel,
    real_t scaleFactor
)
    : FVM::EquationTerm(fluidGrid),
      gridTarget(grid_target),
      gridPrimary(grid_primary),
      unknowns(u),
      id_f_primary(id_f_primary),
      energyKernel(energy_kernel),
      scaleFactor(scaleFactor) {
    SetName("MollerProductionRateTerm");
    Validate();

    id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    AddUnknownForJacobian(unknowns, id_ntot);

    rate.assign(fluidGrid->GetNr(), 0.0);
}

void MollerProductionRateTerm::Validate() const {
    if (grid == nullptr) throw DREAMException("MollerProductionRateTerm: fluid grid is null.");
    if (gridPrimary == nullptr)
        throw DREAMException("MollerProductionRateTerm: primary grid is null.");
    if (gridTarget == nullptr)
        throw DREAMException("MollerProductionRateTerm: target grid is null.");
    if (unknowns == nullptr)
        throw DREAMException("MollerProductionRateTerm: unknown handler is null.");
    if (energyKernel == nullptr)
        throw DREAMException("MollerProductionRateTerm: energyKernel is null.");
    if (grid->GetNr() != gridPrimary->GetNr())
        throw DREAMException(
            "MollerProductionRateTerm: fluid grid and primary grid must have same Nr."
        );
}

void MollerProductionRateTerm::Rebuild(
    const real_t t, const real_t /*dt*/, FVM::UnknownQuantityHandler * /*uqh*/
) {
    // explicit-time cache: rebuild only when time changes beyond roundoff
    constexpr real_t TIME_EPS_FACTOR = 100;
    if (std::abs(t - t_rebuilt) <= TIME_EPS_FACTOR * std::numeric_limits<real_t>::epsilon()) {
        return;
    }
    const real_t *fP = unknowns->GetUnknownDataPrevious(id_f_primary);
    if (fP == nullptr) {
        throw DREAMException("MollerProductionRateTerm: f_primary(previous) is null.");
    }

    len_t offP = 0;
    const len_t Nr = gridPrimary->GetNr();

    for (len_t ir = 0; ir < Nr; ++ir) {
        const auto *mgP = gridPrimary->GetMomentumGrid(ir);
        const len_t Np1P = mgP->GetNp1();
        const len_t NxiP = mgP->GetNp2();

        const real_t *VpP = gridPrimary->GetVp(ir);
        const real_t VpVol = gridPrimary->GetVpVol(ir);
        real_t sum = 0.0;
        for (len_t l = 0; l < NxiP; ++l) {
            const real_t dxi1 = mgP->GetDp2(l);
            const len_t base = offP + l * Np1P;
            for (len_t k = 0; k < Np1P; ++k) {
                const real_t dp1 = mgP->GetDp1(k);
                const len_t idx = base + k;

                const real_t f = fP[idx];
                const real_t Vp1 = VpP[l * Np1P + k];

                // TotalCS is v1*sigma_tot(p1)
                const real_t sigmaTotV = energyKernel->TotalCS(k);

                sum += dp1 * dxi1 * Vp1 / VpVol * f * sigmaTotV;
            }
        }

        rate[ir] = scaleFactor * sum;
        offP += mgP->GetNCells();
    }
    t_rebuilt = t;
}

void MollerProductionRateTerm::SetVectorElements(real_t *vec, const real_t *n_tot) {
    for (len_t ir = 0; ir < grid->GetNr(); ++ir) {
        vec[ir] += n_tot[ir] * rate[ir];
    }
}

bool MollerProductionRateTerm::SetJacobianBlock(
    const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t * /*x*/
) {
    if (derivId != id_ntot) {
        return false;
    }
    // Intended usage: equation row is fluid (n_re or S_particle), derivative w.r.t n_tot.
    // The coupling is local in radius: d/dn_tot[ir] (n_tot[ir]*rate[ir]) = rate[ir].
    for (len_t ir = 0; ir < grid->GetNr(); ++ir) {
        jac->SetElement(ir, ir, rate[ir]);
    }

    return true;
}
