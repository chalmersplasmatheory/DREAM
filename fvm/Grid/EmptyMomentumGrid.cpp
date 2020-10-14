/**
 * Implementation of an empty momentum grid (to be used
 * with spatial 1D grids).
 */

#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace DREAM::FVM;


/**
 * Generate the empty momentum grid, with np=nxi=1, p=1 and xi=1.
 */
bool EmptyMomentumGridGenerator::Rebuild(
    const real_t /*t*/, const len_t /*ir*/, DREAM::FVM::MomentumGrid *mg,
    const DREAM::FVM::RadialGrid* /*rg*/
) {
    const len_t N = 1;

    real_t
        *p1     = new real_t[N],
        *p1_f   = new real_t[N+1],
        *p2     = new real_t[N],
        *p2_f   = new real_t[N+1],
        *dp1    = new real_t[N],
        *dp2    = new real_t[N],
        *p      = new real_t[N*N],
        *p_f1   = new real_t[(N+1)*N],
        *p_f2   = new real_t[N*(N+1)],
        *gamma      = new real_t[N*N],
        *gamma_f1   = new real_t[(N+1)*N],
        *gamma_f2   = new real_t[N*(N+1)],
        *xi0    = new real_t[N*N],
        *xi0_f1 = new real_t[(N+1)*N],
        *xi0_f2 = new real_t[N*(N+1)];

    p1[0] = p1_f[0] = p1_f[1] = dp1[0] = dp2[0] = 1;
    p2_f[0] = p2_f[1] = p2[0] = 1;
    mg->InitializeP1("p",  N, p1, p1_f, dp1, nullptr);
    mg->InitializeP2("xi", N, p2, p2_f, dp2, nullptr);

    // Initialize p and xi0 grids
    for (len_t j = 0; j < N; j++) {
        for (len_t i = 0; i < N; i++) {
            p[j*N + i]   = p1[i];
            gamma[j*N + i] = sqrt(1+p[j*N+i]*p[j*N+i]);
            xi0[j*N + i] = p2[i];
        }
    }

    for (len_t j = 0; j < N; j++) {
        for (len_t i = 0; i < N+1; i++) {
            p_f1[j*(N+1) + i]   = p1_f[i];
            gamma_f1[j*(N+1) + i] = sqrt(1+p_f1[j*(N+1) + i]*p_f1[j*(N+1) + i]);
            xi0_f1[j*(N+1) + i] = p2_f[j];
        }
    }

    for (len_t j = 0; j < N+1; j++) {
        for (len_t i = 0; i < N; i++) {
            p_f2[j*N + i]   = p1[i];
            gamma_f2[j*N + i] = sqrt(1+p_f2[j*N + i]*p_f2[j*N + i]);
            xi0_f2[j*N + i] = p2_f[j];
        }
    }
    
    mg->InitializePAndXi0(p, p_f1, p_f2, gamma,gamma_f1, gamma_f2, xi0, xi0_f1, xi0_f2);

    this->isBuilt = true;
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
   const len_t , const len_t ,
            fluxGridType , 
            const len_t ntheta, const real_t* ,
            const real_t*, real_t *&sqrtg
) const {
    for (len_t it = 0; it < ntheta; it++) {
        sqrtg[it] = 1;
    }
}

