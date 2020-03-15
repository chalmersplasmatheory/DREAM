/**
 * Implementation of a cylindrical radial grid (corresponding to the
 * large-aspect ratio limit for a tokamak).
 */

#include <cmath>
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"


using namespace TQS::FVM;

/**
 * Constructor.
 *
 * nx: Number of radial grid points.
 * B0: Magnetic field strength.
 * x0: Value of inner radial flux grid point.
 * xa: Value of outer radial flux grid point.
 */
CylindricalRadialGridGenerator::CylindricalRadialGridGenerator(
    const len_t nx, const real_t B0,
    const real_t x0, const real_t xa
) : nx(nx), xMin(x0), xMax(xa), B0(B0) {}


/*************************************
 * PUBLIC METHODS                    *
 *************************************/
/**
 * (Re-)builds the given radial grid.
 *
 * rGrid: Radial grid to re-build.
 */
bool CylindricalRadialGridGenerator::Rebuild(const real_t, RadialGrid *rGrid) {
    real_t
        *x    = new real_t[nx],
        *x_f  = new real_t[nx+1],
        *dx   = new real_t[nx],
        *dx_f = new real_t[nx-1];

    real_t
        *volumes       = new real_t[nx],
        *avGradr2      = new real_t[nx],
        *avGradr2_R2_f = new real_t[nx+1];

    // Construct flux grid
    for (len_t i = 0; i < nx; i++)
        dx[i] = (xMax - xMin) / nx;

    for (len_t i = 0; i < nx+1; i++)
        x_f[i] = xMin + i*dx[0];

    // Construct cell grid
    for (len_t i = 0; i < nx; i++)
        x[i] = 0.5 * (x_f[i+1] + x_f[i]);

    for (len_t i = 0; i < nx-1; i++)
        dx_f[i] = x[i+1] - x[i];

    // Construct grid volumes
    for (len_t i = 0; i < nx; i++)
        volumes[i] = 4*M_PI*x[i]*dx[i];

    // Construct jacobians
    for (len_t i = 0; i < nx; i++)
        avGradr2[i] = x[i];

    for (len_t i = 0; i < nx+1; i++)
        avGradr2_R2_f[i] = x_f[i];

    rGrid->Initialize(
        x, x_f, dx, dx_f,
        volumes, avGradr2, avGradr2_R2_f
    );

    this->isBuilt = true;

    return true;
}

