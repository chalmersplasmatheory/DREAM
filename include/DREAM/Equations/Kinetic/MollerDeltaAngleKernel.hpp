#ifndef _DREAM_EQUATIONS_MOLLER_DELTA_ANGLE_KERNEL_HPP
#define _DREAM_EQUATIONS_MOLLER_DELTA_ANGLE_KERNEL_HPP

#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/config.h"

namespace DREAM {

/**
 * Helper representing the orbit-averaged kinematic delta function
 * appearing in the Møller knock-on problem.
 */
class MollerDeltaAngleKernel {
    struct XiStarInterp {
        len_t m0;
        real_t w1;
    };

    struct DeltaInterp {
        const real_t *D0;  // slice 0 (NxiK x NxiP), column-major in l
        const real_t *D1;  // slice 1 or == D0 if clamped
        real_t w0, w1;     // interpolation weights (w0+w1=1), (1,0) if clamped
    };

   private:
    const FVM::Grid *gridK = nullptr;
    const FVM::Grid *gridP = nullptr;

    real_t pCutoff = 0;
    len_t nXiStarsTabulate = 0;
    len_t nPointsIntegral = 0;
    KnockOnUtilities::orbit_integration_method integrationMethod;

    // xi_star tabulation grid
    real_t xiStarMin = 0, xiStarMax = 1, dXiStar = 1;
    real_t *xiStarsTab = nullptr;

    // deltaTable[ir][m] -> slice pointer, each plane is contiguous storage of size NxiK*NxiP
    real_t ***deltaTable = nullptr;
    real_t **deltaTableStorage = nullptr;

    // interpolation metadata per (ir,i,k)
    XiStarInterp **xiInterp = nullptr;

    void Deallocate();
    void ValidateGridAssumptions() const;
    void ValidateInputParameters() const;

    void BuildXiStarTabulationGrid();
    void AllocateDeltaTables();
    void TabulateDeltaMatrixOnXiStarGrid();
    void EvaluateXiStarRangeOnGrids(real_t &xiStarMin, real_t &xiStarMax);
    void BuildXiStarInterp();

    const DeltaInterp GetDeltaInterp(len_t ir, len_t i, len_t k) const;

   public:
    MollerDeltaAngleKernel(
        const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t p_cutoff,
        len_t n_xi_stars_tabulate, len_t n_points_integral,
        KnockOnUtilities::orbit_integration_method integration_method
    );
    ~MollerDeltaAngleKernel();

    void GridRebuilt();

    void AccumulatePitch(
        len_t ir, len_t i, len_t k, const real_t *W_l, real_t Sik, real_t *outPitch_j
    ) const;
};

}  // namespace DREAM

#endif
