#ifndef _DREAM_EQUATIONS_MOLLER_DELTA_ANGLE_KERNEL_HPP
#define _DREAM_EQUATIONS_MOLLER_DELTA_ANGLE_KERNEL_HPP

#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "DREAM/Equations/Kinetic/MollerPitchContraction.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/config.h"

namespace DREAM {

/**
 * Normalized angular shape for Møller knock-ons, averaged over both
 * outgoing and primary pitch cells with the continuum orbit measure.
 */
class MollerDeltaAngleKernel {
    using XiStarInterp = MollerPitchInterpolation;

   private:
    const FVM::Grid *gridK = nullptr;
    const FVM::Grid *gridP = nullptr;

    real_t pCutoff = 0;
    len_t nXiStarsTabulate = 0;
    len_t poloidalOrder = 0, pitchOrder = 0;
    OptionConstants::boltzmann_model model;


    // xi_star tabulation grid
    real_t xiStarMin = 0, xiStarMax = 1, dXiStar = 1;
    real_t *xiStarsTab = nullptr;
    bool hasActivePairs = false;

    // Each slice stores NxiK entries per primary column.
    real_t ***deltaTable = nullptr;
    real_t **deltaTableStorage = nullptr;
    len_t nrAllocated = 0;

    // interpolation metadata per (ir,i,k)
    XiStarInterp **xiInterp = nullptr;

    void Deallocate();
    void ValidateGridAssumptions() const;
    void ValidateInputParameters() const;

    void BuildXiStarTabulationGrid();
    void AllocateDeltaTables();
    void TabulateDeltaMatrixOnXiStarGrid();
    void EvaluateXiStarRangeOnGrids(real_t &xiStarMin, real_t &xiStarMax);
    bool EvaluateXiStarOnGrids(len_t i, len_t k, real_t &xiStar) const;
    void BuildXiStarInterp();


   public:
    MollerDeltaAngleKernel(
        const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t p_cutoff,
        len_t n_xi_stars_tabulate, len_t poloidal_order, len_t pitch_order=4,
        OptionConstants::boltzmann_model model=OptionConstants::BOLTZMANN_MODEL_FULL_MOLLER
    );
    ~MollerDeltaAngleKernel();

    void GridRebuilt();

    bool UsesChiuHarvey() const { return model == OptionConstants::BOLTZMANN_MODEL_CHIU_HARVEY; }
    len_t GetPrimaryPitchColumns(len_t ir) const { return UsesChiuHarvey() ? 2 : gridP->GetNp2(ir); }
    len_t GetPitchScratchSize() const { return nXiStarsTabulate * GetPrimaryPitchColumns(0); }

    // In Chiu–Harvey mode, each momentum has weights for the +1 and -1 poles, in that order.
    void AccumulatePitch(
        len_t ir, len_t i, len_t k, const real_t *W_l, real_t Sik, real_t *outPitch_j
    ) const;

    // Batched interface used by the kinetic source; scratch is caller-owned.
    void AccumulatePitchAllPrimaries(
        len_t ir, len_t i, const real_t *W_kl, const real_t *S_ik,
        real_t *scratch, real_t *outPitch_j
    ) const;
};

}  // namespace DREAM

#endif
