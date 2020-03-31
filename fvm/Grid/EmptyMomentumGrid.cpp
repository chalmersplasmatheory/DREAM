/**
 * Implementation of an empty momentum grid (to be used
 * with spatial 1D grids).
 */

#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace DREAM::FVM;


/**
 * Generate the empty momentum grid.
 */
bool EmptyMomentumGridGenerator::Rebuild(
    const real_t t, const len_t ir, DREAM::FVM::MomentumGrid *mg,
    const DREAM::FVM::RadialGrid *rg
) {
    const len_t N = 1;

    real_t
        *p1   = new real_t[N],
        *p1_f = new real_t[N+1],
        *p2   = new real_t[N],
        *p2_f = new real_t[N+1],
        *dp1  = new real_t[N],
        *dp2  = new real_t[N];

    p1[0] = p2[0] = dp1[0] = dp2[0] = 0;
    p1_f[0] = p1_f[1] = p2_f[0] = p2_f[1] = 0;

    mg->InitializeP1(N, p1, p1_f, dp1, nullptr);
    mg->InitializeP2(N, p2, p2_f, dp2, nullptr);
    
    return true;
}

/**
 * Evaluate the metric for this grid.
 *
 * p1, p2:    [Unused]
 * ir:        Index of radial grid point to evaluate metric in.
 * rGrid:     Radial grid used for evaluating metric.
 * ntheta:    Number of poloidal angle grid points.
 * theta:     Poloidal angle grid.
 * rFluxGrid: If 'true', indicates that the calculation
 *            should be done on the radial flux grid.
 *
 * RETURNS
 * sqrtg: Square root of the metric trace.
 */
void EmptyMomentumGrid::EvaluateMetric(
    const real_t, const real_t,
    const len_t ir, const RadialGrid *rGrid,
    const len_t ntheta, const real_t *theta,
    bool rFluxGrid, real_t *sqrtg
) const {
    const real_t *B = (
        rFluxGrid ?
            rGrid->BOfTheta_f(ir) :
            rGrid->BOfTheta(ir)
    );

    // TODO calculate correctly!
    for (len_t i = 0; i < ntheta; i++) {
        sqrtg[i] = 1;
    }
}

