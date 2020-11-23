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
) : RadialGridGenerator(nx), xMin(x0), xMax(xa), B0(B0) {
    isUpDownSymmetric = true;
    ntheta_interp = 1;
}


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
    x    = new real_t[GetNr()];   
    x_f  = new real_t[GetNr()+1];
    real_t
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

    // Construct magnetic field quantities
    R0            = std::numeric_limits<real_t>::infinity();
    BtorGOverR0   = new real_t[GetNr()];
    psiPrimeRef   = new real_t[GetNr()];
    BtorGOverR0_f = new real_t[GetNr()+1];
    psiPrimeRef_f = new real_t[GetNr()+1];    
    for (len_t ir = 0; ir < GetNr(); ir++){
        BtorGOverR0[ir] = B0;
        psiPrimeRef[ir] = 0; // no poloidal magnetic field; the result of including would only be a radially dependent constant added to B0
    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        BtorGOverR0_f[ir] = B0;
        psiPrimeRef_f[ir] = 0; // no poloidal magnetic field
    }
    rGrid->SetReferenceMagneticFieldData(
        BtorGOverR0, BtorGOverR0_f, psiPrimeRef, psiPrimeRef_f, R0
    );

    this->isBuilt = true;

    return true;
}
