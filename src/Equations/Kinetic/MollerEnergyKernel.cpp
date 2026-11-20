#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"

using namespace DREAM;

MollerEnergyKernel::MollerEnergyKernel(
    const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t p_cutoff
)
    : gridK(grid_knockon), gridP(grid_primary), pCutoff(p_cutoff) {
    ValidateInputParameters();
    GridRebuilt();
}

MollerEnergyKernel::~MollerEnergyKernel() { Deallocate(); }

void MollerEnergyKernel::Deallocate() {
    if (Sik != nullptr) {
        delete[] Sik;
        Sik = nullptr;
    }
    NpK0 = 0;
    NpP0 = 0;
}

void MollerEnergyKernel::ValidateInputParameters() const {
    if (!(pCutoff > 0)) {
        throw DREAMException("MollerEnergyKernel: invalid pCutoff=%.16g (must be > 0).", pCutoff);
    }
}

void DREAM::MollerEnergyKernel::ValidateGridAssumptions() const {
    if (gridK == nullptr || gridP == nullptr)
        throw DREAMException("MollerEnergyKernel: grid pointers must not be null.");

    if (!(pCutoff > 0))
        throw DREAMException("MollerEnergyKernel: invalid pCutoff=%.16g (must be > 0).", pCutoff);

    const len_t NrK = gridK->GetNr();
    const len_t NrP = gridP->GetNr();
    if (NrK != NrP) throw DREAMException("MollerEnergyKernel: gridK and gridP must have same Nr.");

    // Enforce uniform momentum resolution across radii (same as operator assumption).
    for (len_t ir = 1; ir < NrK; ++ir) {
        if (gridK->GetNp1(ir) != gridK->GetNp1(0) || gridK->GetNp2(ir) != gridK->GetNp2(0)) {
            throw DREAMException(
                "MollerEnergyKernel: requires uniform knock-on momentum grid across radii."
            );
        }
        if (gridP->GetNp1(ir) != gridP->GetNp1(0) || gridP->GetNp2(ir) != gridP->GetNp2(0)) {
            throw DREAMException(
                "MollerEnergyKernel: requires uniform primary momentum grid across radii."
            );
        }
    }
}

void MollerEnergyKernel::GridRebuilt() {
    Deallocate();
    ValidateGridAssumptions();

    NpK0 = gridK->GetNp1(0);
    NpP0 = gridP->GetNp1(0);

    Sik = new real_t[NpK0 * NpP0];

    for (len_t i = 0; i < NpK0; ++i) {
        for (len_t k = 0; k < NpP0; ++k) {
            Sik[i * NpP0 + k] = KnockOnUtilities::EvaluateMollerFluxMatrixElementOnGrid(
                i, k, gridK, gridP, pCutoff
            );
        }
    }
}

// Velocity * total cross section (v_1 * sigma_tot(p_1)) at incident momentum p1_k
real_t MollerEnergyKernel::TotalCS(len_t k) const {
    return KnockOnUtilities::EvaluateMollerFluxIntegratedOverKnockonGrid(k, gridK, gridP, pCutoff);
}
