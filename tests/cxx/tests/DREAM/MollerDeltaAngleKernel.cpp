#include "MollerDeltaAngleKernel.hpp"

#include <cmath>
#include <vector>

#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAMTESTS::_DREAM;

bool MollerDeltaAngleKernel::Run(bool) {
    bool success = true;

    if (CheckMDK_DeltaInterpolationConservation())
        this->PrintOK("Interpolated delta columns preserve normalization.");
    else {
        success = false;
        this->PrintError("Test failed: interpolated delta columns do not normalize to 1.");
    }
    return success;
}

bool MollerDeltaAngleKernel::CheckMDK_DeltaInterpolationConservation() {
    const real_t tol = 1e-14;

    len_t nr = 2;
    len_t npK = 5;
    len_t nxiK = 40;
    len_t npP = 7;
    len_t nxiP = 30;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, npK, nxiK, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, npP, nxiP, ntheta_interp, nrProfiles, pMin, pMax);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1(1);
    constexpr len_t n_xi_stars_tabulate = 80;
    constexpr len_t n_points_integral = 80;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    DREAM::MollerDeltaAngleKernel K(
        gridK, gridP, pCutoff, n_xi_stars_tabulate, n_points_integral, integrationMethod
    );
    K.GridRebuilt();

    bool success = true;

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mgK = gridK->GetMomentumGrid(ir);
        auto *mgP = gridP->GetMomentumGrid(ir);

        const len_t NpK = mgK->GetNp1();
        const len_t NxiK = mgK->GetNp2();
        const len_t NpP = mgP->GetNp1();
        const len_t NxiP = mgP->GetNp2();

        std::vector<real_t> W_l(NxiP, 0.0);
        std::vector<real_t> outPitch(NxiK, 0.0);

        for (len_t i = 0; i < NpK; i++) {
            for (len_t k = 0; k < NpP; k++) {
                for (len_t l = 0; l < NxiP; l++) {
                    std::fill(W_l.begin(), W_l.end(), 0.0);
                    W_l[l] = 1.0;

                    std::fill(outPitch.begin(), outPitch.end(), 0.0);
                    K.AccumulatePitch(ir, i, k, W_l.data(), /*Sik=*/1.0, outPitch.data());

                    real_t sum = 0.0;
                    for (len_t j = 0; j < NxiK; j++) sum += mgK->GetDp2(j) * outPitch[j];

                    const bool primaryVoid = (gridP->GetVpOverP2AtZero(ir)[l] == 0);
                    if (primaryVoid) {
                        if (fabs(sum) != 0) {
                            this->PrintError(
                                "Delta normalization non-zero in void primary cell at ir=%ld i=%ld "
                                "k=%ld l=%ld: sum=%.16g",
                                ir, i, k, l, sum
                            );
                            success = false;
                        }
                    } else {
                        if (fabs(sum - 1.0) > tol) {
                            this->PrintError(
                                "Delta normalization failed at ir=%ld i=%ld k=%ld l=%ld: sum=%.16g",
                                ir, i, k, l, sum
                            );
                            success = false;
                        }
                    }
                }
            }
        }
    }

    delete gridF;
    delete gridK;
    delete gridP;

    return success;
}
