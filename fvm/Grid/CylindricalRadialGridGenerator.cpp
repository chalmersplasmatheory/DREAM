/**
 * Implementation of a cylindrical radial grid (corresponding to the
 * large-aspect ratio limit for a tokamak).
 */

#include <cmath>
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"


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
    const len_t nx, const real_t B0,
    const real_t x0, const real_t xa
) : RadialGridGenerator(nx), xMin(x0), xMax(xa), B0(B0) {}


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

    // Construct magnetic field quantities
    len_t ntheta = 1;
    real_t
        *theta = new real_t[ntheta],
        *B     = new real_t[GetNr()*ntheta],
        *B_f   = new real_t[(GetNr()+1)*ntheta];

    theta[0] = 0;

    for (len_t i = 0; i < GetNr(); i++)
        B[i] = B0;
    for (len_t i = 0; i < GetNr()+1; i++)
        B_f[i] = B0;

    rGrid->InitializeMagneticField(
        ntheta, theta, B, B_f
    );

    this->isBuilt = true;

    return true;
}

/**
 * Re-build the phase space jacobians.
 *
 * rGrid: Radial grid to build jacobians for.
 */
void CylindricalRadialGridGenerator::RebuildJacobians(RadialGrid *rGrid) {
    real_t
        **Vp      = new real_t*[GetNr()],
        **Vp_fr   = new real_t*[GetNr()+1],
        **Vp_f1   = new real_t*[GetNr()],
        **Vp_f2   = new real_t*[GetNr()];

    // Poloidal angle grid
    real_t theta = 0;

    // Set Vp, Vp_f1 and Vp_f2
    for (len_t ir = 0; ir < GetNr(); ir++) {
        const MomentumGrid *mg = rGrid->GetMomentumGrid(ir);
        const real_t r = rGrid->GetR(ir);
        const real_t J = 4*M_PI*r;

        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();
        const real_t
            *p1   = mg->GetP1(),
            *p2   = mg->GetP2(),
            *p1_f = mg->GetP1_f(),
            *p2_f = mg->GetP2_f();

        Vp[ir]    = new real_t[n1*n2];
        Vp_f1[ir] = new real_t[(n1+1)*n2];
        Vp_f2[ir] = new real_t[n1*(n2+1)];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                real_t v;
                mg->EvaluateMetric(p1[i], p2[j], ir, rGrid, 1, &theta, false, &v);

                Vp[ir][j*n1 + i] = J*v;
            }
        }

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1+1; i++) {
                real_t v;
                mg->EvaluateMetric(p1_f[i], p2[j], ir, rGrid, 1, &theta, false, &v);

                Vp_f1[ir][j*(n1+1) + i] = J*v;
            }
        }

        for (len_t j = 0; j < n2+1; j++) {
            for (len_t i = 0; i < n1; i++) {
                real_t v;
                mg->EvaluateMetric(p1[i], p2_f[j], ir, rGrid, 1, &theta, false, &v);

                Vp_f2[ir][j*n1 + i] = J*v;
            }
        }
    }

    // Set Vp_fr
    for (len_t ir = 0; ir < GetNr()+1; ir++) {
        // We inherently assume that the momentum grids at all
        // radii are the same here, so we might as well just evaluate
        // everything on the innermost momentum grid
        const MomentumGrid *mg = rGrid->GetMomentumGrid(0);
        const real_t r = rGrid->GetR_f(ir);
        const real_t J = 4*M_PI*r;

        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();
        const real_t
            *p1   = mg->GetP1(),
            *p2   = mg->GetP2(),
            *p1_f = mg->GetP1_f(),
            *p2_f = mg->GetP2_f();

        Vp_fr[ir] = new real_t[n1*n2];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                real_t v;
                mg->EvaluateMetric(p1[i], p2[j], ir, rGrid, 1, &theta, true, &v);

                Vp_fr[ir][j*n1 + i] = J*v;
            }
        }
    }

    rGrid->InitializeVprime(Vp, Vp_fr, Vp_f1, Vp_f2);
}

