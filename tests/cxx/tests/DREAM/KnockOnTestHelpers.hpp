#ifndef _DREAMTESTS_DREAM_KNOCK_ON_TEST_HELPERS_HPP
#define _DREAMTESTS_DREAM_KNOCK_ON_TEST_HELPERS_HPP

#include <vector>

#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAMTESTS::_DREAM {

inline void SetPreviousUnknownData(
    DREAM::FVM::UnknownQuantityHandler *uqh, const len_t id, const real_t *data,
    const real_t t = 0.0
) {
    uqh->Store(id, data, /*offs*/ 0, /*mayBeConstant*/ false);
    uqh->SaveStep(t, /*trueSave*/ false);
}

inline DREAM::FVM::UnknownQuantityHandler *BuildUQH_Minimal(
    DREAM::FVM::Grid *grid_fluid, DREAM::FVM::Grid *grid_primary, len_t &id_ntot, len_t &id_E,
    len_t &id_f_primary
) {
    auto *uqh = new DREAM::FVM::UnknownQuantityHandler();

    id_ntot = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, "0", grid_fluid);
    id_E = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, "0", grid_fluid);
    id_f_primary = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_HOT, "0", grid_primary);

    // Initial values for n_tot and E_field
    real_t *tmp = new real_t[grid_fluid->GetNr()];
    for (len_t ir = 0; ir < grid_fluid->GetNr(); ir++) tmp[ir] = 1.0;
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_N_TOT, tmp);

    for (len_t ir = 0; ir < grid_fluid->GetNr(); ir++) tmp[ir] = 0.0;
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_E_FIELD, tmp);
    delete[] tmp;

    // f_primary initial value (zeros)
    {
        const len_t N = grid_primary->GetNCells();
        real_t *f0 = new real_t[N];
        for (len_t idx = 0; idx < N; idx++) f0[idx] = 0.0;
        uqh->SetInitialValue(DREAM::OptionConstants::UQTY_F_HOT, f0);
        delete[] f0;
    }

    return uqh;
}

inline DREAM::FVM::UnknownQuantityHandler *BuildUQH_MinimalFRe(
    DREAM::FVM::Grid *grid_fluid, DREAM::FVM::Grid *grid_runaway, len_t &id_ntot, len_t &id_E,
    len_t &id_f_re
) {
    auto *uqh = new DREAM::FVM::UnknownQuantityHandler();

    id_ntot = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, "0", grid_fluid);
    id_E = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, "0", grid_fluid);
    id_f_re = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_RE, "0", grid_runaway);

    real_t *tmp = new real_t[grid_fluid->GetNr()];
    for (len_t ir = 0; ir < grid_fluid->GetNr(); ir++) tmp[ir] = 1.0;
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_N_TOT, tmp);

    for (len_t ir = 0; ir < grid_fluid->GetNr(); ir++) tmp[ir] = 0.0;
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_E_FIELD, tmp);
    delete[] tmp;

    {
        const len_t N = grid_runaway->GetNCells();
        real_t *f0 = new real_t[N];
        for (len_t idx = 0; idx < N; idx++) f0[idx] = 0.0;
        uqh->SetInitialValue(DREAM::OptionConstants::UQTY_F_RE, f0);
        delete[] f0;
    }

    return uqh;
}

inline real_t IntegrateTotalProductionOverKnockonGrid(
    const DREAM::FVM::Grid *grid_knockon, const real_t *sourceVector
) {
    real_t total = 0.0;
    len_t offset = 0;
    for (len_t ir = 0; ir < grid_knockon->GetNr(); ir++) {
        auto *mg = grid_knockon->GetMomentumGrid(ir);
        const len_t Np = mg->GetNp1();
        const len_t Nxi = mg->GetNp2();

        for (len_t i = 0; i < Np; i++) {
            const real_t dp = mg->GetDp1(i);
            for (len_t j = 0; j < Nxi; j++) {
                const real_t dxi = mg->GetDp2(j);
                const len_t ind = offset + Np * j + i;
                const real_t Vp = grid_knockon->GetVp(ir, i, j);
                total += dp * dxi * Vp * sourceVector[ind];
            }
        }
        offset += mg->GetNCells();
    }
    return total;
}

inline real_t PredictTotalProductionFromMollerS(
    const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary,
    const real_t *f_primary, real_t pCutoff
) {
    real_t total = 0.0;
    len_t offsetP = 0;

    for (len_t ir = 0; ir < grid_primary->GetNr(); ir++) {
        auto *mgP = grid_primary->GetMomentumGrid(ir);

        for (len_t k = 0; k < mgP->GetNp1(); k++) {
            const real_t dp1 = mgP->GetDp1(k);

            const real_t sigmaTot =
                DREAM::KnockOnUtilities::EvaluateMollerFluxIntegratedOverKnockonGrid(
                    k, grid_knockon, grid_primary, pCutoff
                );

            for (len_t l = 0; l < mgP->GetNp2(); l++) {
                const real_t dxi1 = mgP->GetDp2(l);
                const len_t indP = offsetP + mgP->GetNp1() * l + k;
                const real_t f = f_primary[indP];
                const real_t Vp1 = grid_primary->GetVp(ir, k, l);
                total += dp1 * dxi1 * Vp1 * f * sigmaTot;
            }
        }

        offsetP += mgP->GetNCells();
    }

    return total;
}

}  // namespace DREAMTESTS::_DREAM

#endif
