/**
 * Implementation of a cylindrical radial grid (corresponding to the
 * large-aspect ratio limit for a tokamak).
 */

#include <cmath>
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include <functional>

using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * nx: Number of radial grid points.
 * B0: Magnetic field strength.
 * x0: Value of inner radial flux grid point.
 * xa: Value of outer radial flux grid point.
 */
CylindricalRadialGridGenerator::CylindricalRadialGridGenerator(
     len_t nx,  real_t B0,
     real_t x0, real_t xa
) : RadialGridGenerator(nx), xMin(x0), xMax(xa), B0(B0) {}


/*************************************
 * PUBLIC METHODS                    *
 *************************************/
/**
 * (Re-)builds the given radial grid and creates
 * magnetic field and (spatial) Jacobian data.
 *
 * rGrid: Radial grid to re-build.
 */
bool CylindricalRadialGridGenerator::Rebuild(const real_t, RadialGrid *rGrid) {
    real_t
        *x    = new real_t[GetNr()],
        *x_f  = new real_t[GetNr()+1],
        *dx   = new real_t[GetNr()],
        *dx_f = new real_t[GetNr()-1];

    // Construct flux grid
    for (len_t i = 0; i < GetNr(); i++)
        dx[i] = (xMax - xMin) / GetNr();

    for (len_t i = 0; i < GetNr()+1; i++)
        x_f[i] = xMin + i*dx[0];

    // Construct cell grid
    for (len_t i = 0; i < GetNr(); i++)
        x[i] = 0.5 * (x_f[i+1] + x_f[i]);

    for (len_t i = 0; i < GetNr()-1; i++)
        dx_f[i] = x[i+1] - x[i];

    rGrid->Initialize(x, x_f, dx, dx_f);

    

    this->isBuilt = true;

    return true;
}

/**
 * Re-build the phase space jacobians.
 *
 * grid: Grid to build jacobians for.
 */
void CylindricalRadialGridGenerator::RebuildJacobians(RadialGrid *rGrid, MomentumGrid **momentumGrids, MagneticQuantityHandler *mgnQtyHandler) {

    CreateMagneticFieldData(rGrid->GetR(),rGrid->GetR_f());

    rGrid->InitializeMagneticField(
        ntheta_ref, theta_ref, B_ref, B_ref_f, Bmin_ref, Bmin_ref_f, Bmax_ref, Bmax_ref_f
    );
    
    mgnQtyHandler->Initialize(momentumGrids,
                     ntheta_ref, theta_ref, 
                     B_ref, Jacobian_ref,
                     ROverR0_ref, NablaR2_ref,
                     B_ref_f, Jacobian_ref_f,
                     ROverR0_ref_f, NablaR2_ref_f);

    rGrid->InitializeVprime(mgnQtyHandler->GetVp(0),mgnQtyHandler->GetVp(1),
                            mgnQtyHandler->GetVp(2),mgnQtyHandler->GetVp(3),
                            mgnQtyHandler->GetVpVol(false), mgnQtyHandler->GetVpVol(true));


    

}

void CylindricalRadialGridGenerator::CreateMagneticFieldData(const real_t *x, const real_t *x_f) {
    // Construct magnetic field quantities
    DeallocateMagneticFieldData();

    ntheta_ref = 1;
    theta_ref = new real_t[ntheta_ref];
    theta_ref[0] = 0;
    B_ref          = new real_t*[GetNr()];
    Jacobian_ref   = new real_t*[GetNr()];
    ROverR0_ref    = new real_t*[GetNr()];
    NablaR2_ref    = new real_t*[GetNr()];
    Bmin_ref       = new real_t[GetNr()];
    Bmax_ref       = new real_t[GetNr()];
    B_ref_f        = new real_t*[(GetNr()+1)];
    Jacobian_ref_f = new real_t*[(GetNr()+1)];
    ROverR0_ref_f  = new real_t*[(GetNr()+1)];
    NablaR2_ref_f  = new real_t*[(GetNr()+1)];
    Bmin_ref_f     = new real_t[GetNr()+1];
    Bmax_ref_f     = new real_t[GetNr()+1];
    

    for (len_t ir = 0; ir < GetNr(); ir++){
        B_ref[ir] = new real_t[ntheta_ref];
        Jacobian_ref[ir] = new real_t[ntheta_ref];
        ROverR0_ref[ir]  = new real_t[ntheta_ref];
        NablaR2_ref[ir]  = new real_t[ntheta_ref];
        
        B_ref[ir][0] = B0;
        Bmin_ref[ir] = B0;
        Bmax_ref[ir] = B0;
        Jacobian_ref[ir][0] = x[ir];
        ROverR0_ref[ir][0]  = 1;
        NablaR2_ref[ir][0]  = 1;
    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        B_ref_f[ir] = new real_t[ntheta_ref];
        Jacobian_ref_f[ir] = new real_t[ntheta_ref];
        ROverR0_ref_f[ir]  = new real_t[ntheta_ref];
        NablaR2_ref_f[ir]  = new real_t[ntheta_ref];

        B_ref_f[ir][0] = B0;
        Bmin_ref_f[ir] = B0;
        Bmax_ref_f[ir] = B0;
        Jacobian_ref_f[ir][0] = x_f[ir];
        ROverR0_ref_f[ir][0]  = 1;
        NablaR2_ref_f[ir][0]  = 1;
    }

    
}


void CylindricalRadialGridGenerator::DeallocateMagneticFieldData(){
    if (B_ref==nullptr)
        return;

    for(len_t ir = 0; ir<GetNr(); ir++){
        delete [] B_ref[ir];
        delete [] Jacobian_ref[ir];
        delete [] ROverR0_ref[ir];
        delete [] NablaR2_ref[ir];
    }
    for(len_t ir = 0; ir<GetNr()+1; ir++){
        delete [] B_ref_f[ir];
        delete [] Jacobian_ref_f[ir];
        delete [] ROverR0_ref_f[ir];
        delete [] NablaR2_ref_f[ir];
    }

    delete [] theta_ref;
    delete [] Bmin_ref;
    delete [] Bmin_ref_f;
    delete [] Bmax_ref;
    delete [] Bmax_ref_f;
    delete [] B_ref;
    delete [] Jacobian_ref;
    delete [] ROverR0_ref;
    delete [] NablaR2_ref;
    delete [] B_ref_f;
    delete [] Jacobian_ref_f;
    delete [] ROverR0_ref_f;
    delete [] NablaR2_ref_f;
}
