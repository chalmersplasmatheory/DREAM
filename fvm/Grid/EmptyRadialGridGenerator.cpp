/**
 * Implementation of a grid with 1 grid point at x=0
 */

#include "FVM/Grid/EmptyRadialGrid.hpp"

using namespace DREAM::FVM;



/*************************************
 * PUBLIC METHODS                    *
 *************************************/
/**
 * (Re-)builds the given radial grid and creates
 * magnetic field and (spatial) Jacobian data.
 *
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

    


    return true;
}


void EmptyRadialGridGenerator::CreateMagneticFieldData(const real_t*, const real_t*) {
    // Construct magnetic field quantities
    
    theta_ref  = new real_t[ntheta_ref];
    R0         = 0;
    B_ref          = new real_t*[GetNr()];
    Jacobian_ref   = new real_t*[GetNr()];
    ROverR0_ref    = new real_t*[GetNr()];
    NablaR2_ref    = new real_t*[GetNr()];
    Bmin           = new real_t[GetNr()];
    Bmax           = new real_t[GetNr()];
    BtorGOverR0    = new real_t[GetNr()];
    B_ref_f        = new real_t*[(GetNr()+1)];
    Jacobian_ref_f = new real_t*[(GetNr()+1)];
    ROverR0_ref_f  = new real_t*[(GetNr()+1)];
    NablaR2_ref_f  = new real_t*[(GetNr()+1)];
    Bmin_f         = new real_t[GetNr()+1];
    Bmax_f         = new real_t[GetNr()+1];
    BtorGOverR0_f  = new real_t[GetNr()+1];
    
    theta_ref[0] = 0;
    theta_ref[1] = 2*M_PI;
    for (len_t ir = 0; ir < GetNr(); ir++){
        B_ref[ir] = new real_t[ntheta_ref];
        Jacobian_ref[ir] = new real_t[ntheta_ref];
        ROverR0_ref[ir]  = new real_t[ntheta_ref];
        NablaR2_ref[ir]  = new real_t[ntheta_ref];

        for(len_t it=0; it<ntheta_ref; it++){        
            B_ref[ir][it] = 0;
            Jacobian_ref[ir][it] = 1;
            ROverR0_ref[ir][it]  = 1;
            NablaR2_ref[ir][it]  = 1;
        }
        Bmin[ir] = 0;
        Bmax[ir] = 0;
        BtorGOverR0[ir] = 0;
    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        B_ref_f[ir] = new real_t[ntheta_ref];
        Jacobian_ref_f[ir] = new real_t[ntheta_ref];
        ROverR0_ref_f[ir]  = new real_t[ntheta_ref];
        NablaR2_ref_f[ir]  = new real_t[ntheta_ref];
        for(len_t it=0; it<ntheta_ref; it++){        
            B_ref_f[ir][it] = 0;
            Jacobian_ref_f[ir][it] = 1;
            ROverR0_ref_f[ir][it]  = 1;
            NablaR2_ref_f[ir][it]  = 1;
        }
        Bmin_f[ir] = 0;
        Bmax_f[ir] = 0;
        BtorGOverR0_f[ir] = 0;
    }

}

