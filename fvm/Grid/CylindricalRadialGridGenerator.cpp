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

/**
 * Constructor.
 *
 * x_f_input: Grid points on the flux grid (e.g. the cell edges)
 * nx: Number of radial grid points (so that nx+1 is the size of x_f_input).
 * B0: Magnetic field strength.
 */
CylindricalRadialGridGenerator::CylindricalRadialGridGenerator(
     const real_t *x_f_input, len_t nx,  real_t B0
) : RadialGridGenerator(nx), xMin(x_f_input[0]), xMax(x_f_input[nx]), B0(B0) {
    isUpDownSymmetric = true;
    ntheta_interp = 1;

    this->xf_provided = new real_t[nx+1];
    for(len_t i=0; i<nx+1; i++)
        this->xf_provided[i] = x_f_input[i];

    delete [] x_f_input;
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

    // if x_f has been provided to constructor, set specified grid, otherwise uniform
    if(xf_provided==nullptr){
        for (len_t i = 0; i < GetNr(); i++)
            dx[i] = (xMax - xMin) / GetNr();
        for (len_t i = 0; i < GetNr()+1; i++)
            x_f[i] = xMin + i*dx[0];
    } else {
        for (len_t i = 0; i < GetNr()+1; i++)
            x_f[i] = xf_provided[i];
        for (len_t i = 0; i < GetNr(); i++)
            dx[i] = x_f[i+1] - x_f[i];
    }
    delete [] xf_provided;
    xf_provided = nullptr; // upon next rebuild, create uniform grid at the new resolution

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

/**
 * Calculate the radial flux label 'r' corresponding to the given
 * point in the cartesian SPI coordinate system (centred on the
 * magnetic axis).
 */
void CylindricalRadialGridGenerator::GetRThetaPhiFromCartesian(
   real_t *r, real_t *theta, real_t *phi, real_t x, real_t y, real_t , real_t, real_t
) {
    *r = sqrt(x*x+y*y);
    *theta = std::atan2(y,x);
    *phi = 0;
}

/**
 * Calculate the gradient of the radial flux label 'r' 
 * in cartesian SPI coordinates.
 */
void CylindricalRadialGridGenerator::GetGradRCartesian(
    real_t *gradRCartesian, real_t , real_t theta, real_t
) {
    gradRCartesian[0]=cos(theta);
    gradRCartesian[1]=sin(theta);
    gradRCartesian[2]=0;
}

/**
 * Calculate the radial flux label 'r' at the point of closest approach
 * to the magnetic axis along straight the line between the points given by
 * (x1, y1, z1) and (x2, y2, z2) in the SPI coordinate system
 */
real_t CylindricalRadialGridGenerator::FindClosestApproach(
    real_t x1, real_t y1, real_t,
    real_t x2, real_t y2, real_t
) {
    real_t tc = -(x1*(x2-x1)+y1*(y2-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    return sqrt((x1+tc*(x2-x1))*(x1+tc*(x2-x1))+(y1+tc*(y2-y1))*(y1+tc*(y2-y1)));
}

