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
 * p:            Momentum to evaluate metric for.
 * xi:           Pitch to evaluate metric for.
 * ir:           Index of radial grid point to evaluate metric in.
 * rGrid:        Grid used for evaluating metric.
 * ntheta:       Number of poloidal angle grid points.
 * theta:        Poloidal angle grid.
 * fluxGridType: 0 = distribution grid, 1 = radial flux grid,
 *               2 = p1 flux grid, 3 = p2 flux grid. 
 *
 * RETURNS
 * sqrtg:  Square root of the metric trace.
 */
void PXiMomentumGrid::EvaluateMetric(
    const len_t i, const len_t j ,
    len_t fluxGridType, 
    const len_t ntheta, const real_t* /*theta*/,
    const real_t* B, real_t Bmin, real_t *sqrtg
) const {

    real_t p,xi0;
    p = GetP1(i);
    xi0 = GetP2(j);
    if (fluxGridType==2) {
        p   = this->GetP1_f(i);
    } else if (fluxGridType==3) {
        xi0 = this->GetP2_f(j);
    }

    // sqrtg defined so that the local number density is n=int(f(p1,p2) sqrt(g) dp1 dp2 )
    real_t sqrtg_const = 2*M_PI*p*p*xi0/Bmin;
    real_t xi2_particle;
    // sqrtg=0 outside of the orbit (for theta outside of the integration domain, 
    // ie where (1-xi^2)B/Bmin > 1)
    for (len_t it = 0; it < ntheta; it++) {
        xi2_particle = 1- (B[it]/Bmin)*(1-xi0*xi0);
        if (xi2_particle < 0)
            sqrtg[it] = 0;
        else  
            sqrtg[it] = sqrtg_const * B[it] / sqrt( xi2_particle );
    }
}

