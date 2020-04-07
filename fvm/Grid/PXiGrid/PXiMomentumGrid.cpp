/**
 * Implementation of routines for the general p/xi momentum grid.
 */

#include <cmath>
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace DREAM::FVM::PXiGrid;


/**
 * Evaluates the phase space metric (sqrt(g))
 * in the specified phase space point, at all
 * given poloidal angles.
 *
 * p:         Momentum to evaluate metric for.
 * xi:        Pitch to evaluate metric for.
 * ir:        Index of radial grid point to evaluate metric in.
 * rGrid:     Grid used for evaluating metric.
 * ntheta:    Number of poloidal angle grid points.
 * theta:     Poloidal angle grid.
 * rFluxGrid: If true, the metric is to be evaluated on the
 *            radial flux grid.
 *
 * RETURNS
 * sqrtg:  Square root of the metric trace.
 */
void PXiMomentumGrid::EvaluateMetric(
    const real_t p, const real_t xi,
    const len_t ir, const RadialGrid *rGrid,
    const len_t ntheta, const real_t* /*theta*/,
    bool rFluxGrid, real_t *sqrtg
) const {
    // Evaluate magnetic field on theta grid
    const real_t *B = (
        rFluxGrid ?
            rGrid->BOfTheta_f(ir) :
            rGrid->BOfTheta(ir)
    );

    real_t lambda = (1-xi*xi) / B[0];
    for (len_t i = 0; i < ntheta; i++) {
        sqrtg[i] = p*p/2 * B[i] / std::sqrt(1 - lambda*B[i]);
    }
}

