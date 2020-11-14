/**
 * Implementation of routines for the general p/xi momentum grid.
 */

#include <cmath>
#include <algorithm>
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
//#include "FVM/Grid/RadialGrid.hpp"


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
    fluxGridType fluxGridType, 
    const len_t ntheta, const real_t* /*theta*/,
    const real_t *BOverBmin, real_t *&sqrtg
) const {

    real_t p,xi0;
    if (fluxGridType==FVM::FLUXGRIDTYPE_P1) {
        p   = this->GetP_f1(i,j);
        xi0 = this->GetXi0_f1(i,j);
    } else if (fluxGridType==FVM::FLUXGRIDTYPE_P2) {
        p   = this->GetP_f2(i,j);
        xi0 = this->GetXi0_f2(i,j);
    } else {
        p   = this->MomentumGrid::GetP(i,j);
        xi0 = this->GetXi0(i,j);
    }

    // sqrtg defined so that the local number density is n=int(f(p1,p2) sqrt(g) dp1 dp2 )
    real_t xiOverXi0;
    for (len_t it = 0; it < ntheta; it++) {
        if(BOverBmin[it]==1){
            xiOverXi0 = 1;
        } else {
            real_t xi0Sq = xi0*xi0;
            xiOverXi0 = sqrt( (1 - BOverBmin[it] * (1-xi0Sq))/xi0Sq );
        }
        sqrtg[it] = 2*M_PI*p*p*BOverBmin[it]/xiOverXi0;
    }
}
