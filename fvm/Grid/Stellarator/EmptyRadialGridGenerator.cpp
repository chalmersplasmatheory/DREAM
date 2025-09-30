/**
 * Implementation of a grid with 1 grid point at x=0
 */

#include "FVM/Grid/EmptyRadialGrid.hpp"

using namespace DREAM::FVM;



/*************************************
 * PUBLIC METHODS                    *
 *************************************/
/**
 * Rebuilds a trivial grid of zeros and size 1. *
 * rGrid: Radial grid to re-build.
 */
bool EmptyRadialGridGenerator::Rebuild(const real_t, RadialGrid *rGrid) {
    len_t N = 1;
    real_t
        *x    = new real_t[N],
        *x_f  = new real_t[N+1],
        *dx   = new real_t[N],
        *dx_f = new real_t[N-1];

    // Construct flux grid
    for (len_t i = 0; i < N; i++)
        dx[i] = 0;

    for (len_t i = 0; i < N+1; i++)
        x_f[i] = 0;

    // Construct cell grid
    for (len_t i = 0; i < N; i++)
        x[i] = 0;

    for (len_t i = 0; i < N-1; i++)
        dx_f[i] = 0;

    rGrid->Initialize(x, x_f, dx, dx_f);

    // Construct trivial B=0 field data
    R0         = std::numeric_limits<real_t>::infinity();
    BtorGOverR0    = new real_t[N];
    psiPrimeRef    = new real_t[N];
    BtorGOverR0_f  = new real_t[N+1];
    psiPrimeRef_f  = new real_t[N+1];    
    for (len_t ir = 0; ir < N; ir++){
        BtorGOverR0[ir] = 0;
        psiPrimeRef[ir] = 0; 
    }
    for (len_t ir = 0; ir < N+1; ir++){
        BtorGOverR0_f[ir] = 0;
        psiPrimeRef_f[ir] = 0;
    }
    rGrid->SetReferenceMagneticFieldData(
        BtorGOverR0, BtorGOverR0_f, psiPrimeRef, psiPrimeRef_f, R0
    );

    isBuilt=true;
    return true;
}
