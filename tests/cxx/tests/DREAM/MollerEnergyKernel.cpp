#include "MollerEnergyKernel.hpp"

#include <cmath>

#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAMTESTS::_DREAM;

bool MollerEnergyKernel::Run(bool) {
    bool success = true;

    if (CheckMEK_MollerKernelConservation())
        this->PrintOK("Moller kernel integrates to total cross section.");
    else {
        success = false;
        this->PrintError("Test failed: Moller kernel does not integrate to total cross section.");
    }

    return success;
}

bool MollerEnergyKernel::CheckMEK_MollerKernelConservation() {
    const real_t tol = 1e-14;

    len_t nr = 2;
    len_t npK = 40;
    len_t nxiK = 9;
    len_t npP = 25;
    len_t nxiP = 11;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, npK, nxiK, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, npP, nxiP, ntheta_interp, nrProfiles, pMin, pMax);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1_f(2);

    DREAM::MollerEnergyKernel K(gridK, gridP, pCutoff);
    K.GridRebuilt();

    bool success = true;

    // Under the operator’s grid assumptions, K is radius-independent (built on ir=0).
    // We still loop ir for symmetry with the old test.
    for (len_t ir = 0; ir < nr; ir++) {
        auto *mgK = gridK->GetMomentumGrid(ir);
        auto *mgP = gridP->GetMomentumGrid(ir);

        const len_t NpK = mgK->GetNp1();
        const len_t NpP = mgP->GetNp1();
        const real_t *dpK = mgK->GetDp1();

        for (len_t k = 0; k < NpP; k++) {
            real_t sum = 0.0;
            for (len_t i = 0; i < NpK; i++) {
                sum += dpK[i] * K.DifferentialCS(i, k);
            }
            // Preferred: compare against kernel’s own TotalCS(k)
            const real_t expectedKernel = K.TotalCS(k);

            // Also compare against legacy reference from utilities (optional but useful)
            const real_t expectedRef =
                DREAM::KnockOnUtilities::EvaluateMollerFluxIntegratedOverKnockonGrid(
                    k, gridK, gridP, pCutoff
                );

            if (fabs(sum - expectedKernel) > tol) {
                this->PrintError("MEK conservation failed (vs TotalCS) at k=%ld", k);
                this->PrintError("sum: %.16g", sum);
                this->PrintError("TotalCS: %.16g", expectedKernel);
                success = false;
            }
            if (fabs(expectedKernel - expectedRef) > tol) {
                this->PrintError("MEK TotalCS mismatch vs utilities at k=%ld", k);
                this->PrintError("TotalCS: %.16g", expectedKernel);
                this->PrintError("ref: %.16g", expectedRef);
                success = false;
            }
        }
    }

    delete gridF;
    delete gridK;
    delete gridP;

    return success;
}
