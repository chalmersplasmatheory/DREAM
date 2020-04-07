/**
 * Implementation of routines for the general p/xi momentum grid.
 */

#include <cmath>
#include <algorithm>
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
    const len_t ntheta, const real_t *theta,
    bool rFluxGrid, real_t *sqrtg
) const {
    // Evaluate magnetic field on theta grid
    const real_t *B = (
        rFluxGrid ?
            rGrid->BOfTheta_f(ir) :
            rGrid->BOfTheta(ir)
    );
    const real_t Bmin = (
        rFluxGrid ?
            rGrid->GetBmin_f(ir) :
            rGrid->GetBmin(ir)
    );
    
    // sqrtg defined so that the local number density is n=int(f(p1,p2) sqrt(g) dp1 dp2 )
    real_t sqrtg_const = 2*M_PI*p*p*xi/Bmin;
    real_t xi2_particle;
    // sqrtg=0 outside of the orbit (for theta outside of the integration domain, 
    // ie where (1-xi^2)B/Bmin > 1)
    for (len_t it = 0; it < ntheta; it++) {
        xi2_particle = 1- (B[it]/Bmin)*(1-xi*xi);
        if (xi2_particle < 0)
            sqrtg[it] = 0;
        else  
        sqrtg[it] = sqrtg_const * B[it] / sqrt( xi2_particle );
    }
}

