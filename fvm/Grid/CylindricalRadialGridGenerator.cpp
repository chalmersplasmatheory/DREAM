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
        ntheta, theta, B, B_f, B, B_f, x, x_f
    );

    this->isBuilt = true;

    return true;
}

/**
 * Re-build the phase space jacobians.
 *
 * grid: Grid to build jacobians for.
 */
void CylindricalRadialGridGenerator::RebuildJacobians(RadialGrid *rGrid, MomentumGrid **momentumGrids) {
    real_t
        **Vp      = new real_t*[GetNr()],
        **Vp_fr   = new real_t*[GetNr()+1],
        **Vp_f1   = new real_t*[GetNr()],
        **Vp_f2   = new real_t*[GetNr()];


    // Poloidal angle grid
    real_t theta = 0;

    // Set Vp, Vp_f1 and Vp_f2
    for (len_t ir = 0; ir < GetNr(); ir++) {
        const MomentumGrid *mg = momentumGrids[ir];
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
        // XXX: We inherently assume that the momentum grids at all
        // radii are the same here, so we might as well just evaluate
        // everything on the innermost momentum grid
        const MomentumGrid *mg = momentumGrids[0];
        const real_t r = rGrid->GetR_f(ir);
        const real_t J = 4*M_PI*r;

        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();
        const real_t
            *p1   = mg->GetP1(),
            *p2   = mg->GetP2();

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

/**
 * Re-build flux surface averaged quantities.
 *
 * grid: Grid to build them on.
 */
void CylindricalRadialGridGenerator::RebuildFSAvgQuantities(RadialGrid *rGrid, MomentumGrid **momentumGrids) {
    real_t
     *effectivePassingFraction  = new real_t[GetNr()],
     *magneticFieldMRS         = new real_t[GetNr()],
     *magneticFieldMRS_f       = new real_t[GetNr()+1],
     *nablaR2OverR2_avg        = new real_t[GetNr()],
     *nablaR2OverR2_avg_f      = new real_t[GetNr()+1],
     *OneOverR2_avg            = new real_t[GetNr()],
     *OneOverR2_avg_f          = new real_t[GetNr()+1],
     **xiBounceAverage_f1      = new real_t*[GetNr()],
     **xiBounceAverage_f2      = new real_t*[GetNr()],
     **xi21MinusXi2OverB2_f1   = new real_t*[GetNr()],
     **xi21MinusXi2OverB2_f2   = new real_t*[GetNr()];

    for (len_t ir = 0; ir < GetNr(); ir++) {
        effectivePassingFraction[ir] = 1;
        magneticFieldMRS[ir]         = B0;
        nablaR2OverR2_avg[ir]        = 1;
        OneOverR2_avg[ir]            = 1;
        const MomentumGrid *mg = momentumGrids[ir];
        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();
        
        xiBounceAverage_f1[ir] = new real_t[(n1+1)*n2];
        xiBounceAverage_f2[ir] = new real_t[n1*(n2+1)];
        xi21MinusXi2OverB2_f1[ir] = new real_t[(n1+1)*n2];
        xi21MinusXi2OverB2_f2[ir] = new real_t[n1*(n2+1)];
         
        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1+1; i++) {
                xiBounceAverage_f1[ir][j*(n1+1)+i]    = BounceAverageQuantity(rGrid, mg, ir, i, j, 2, [](real_t xi, real_t  ){return xi;} );
                xi21MinusXi2OverB2_f1[ir][j*(n1+1)+i] = BounceAverageQuantity(rGrid, mg, ir, i, j, 2, [](real_t xi, real_t BOverBMin ){return xi*xi*(1-xi*xi)/(BOverBMin*BOverBMin);} );
            }
        }
        for (len_t j = 0; j < n2+1; j++) {
            for (len_t i = 0; i < n1; i++) {
                xiBounceAverage_f2[ir][j*n1+i]    = BounceAverageQuantity(rGrid, mg, ir, i, j, 3, [](real_t xi, real_t  ){return xi;} );
                xi21MinusXi2OverB2_f2[ir][j*n1+i] = BounceAverageQuantity(rGrid, mg, ir, i, j, 3, [](real_t xi, real_t BOverBMin ){return xi*xi*(1-xi*xi)/(BOverBMin*BOverBMin);} );
            }
        }


    }


    for (len_t ir = 0; ir < GetNr()+1; ir++) {
        magneticFieldMRS_f[ir]  = B0;
        nablaR2OverR2_avg_f[ir] = 1;
        OneOverR2_avg_f[ir]     = 1;
    }
    rGrid->InitializeFSAvg(effectivePassingFraction, magneticFieldMRS,magneticFieldMRS_f,
                    xiBounceAverage_f1, xiBounceAverage_f2,
                    xi21MinusXi2OverB2_f1, xi21MinusXi2OverB2_f2,
                    nablaR2OverR2_avg, nablaR2OverR2_avg_f,
                    OneOverR2_avg, OneOverR2_avg_f);
}




real_t CylindricalRadialGridGenerator::BounceAverageQuantity(RadialGrid *, const MomentumGrid *mg,  len_t , len_t i, len_t j, len_t fluxGrid , std::function<real_t(real_t,real_t)> F){
    real_t xi0;
    if (fluxGrid==2)
        xi0 = mg->GetXi0_f1(i,j);
    else if (fluxGrid==3)
        xi0 = mg->GetXi0_f2(i,j);
    else 
        xi0 = mg->GetXi0(i,j);
    return F(xi0,1);
}
real_t CylindricalRadialGridGenerator::FluxSurfaceAverageQuantity(RadialGrid *,len_t , bool , std::function<real_t(real_t)> F){
    return F(1);
}


/** 
 * Sketching the flux surface averaging function of an arbitrary quantity F = F(xi,B) 

// Calculates the bounce average of an arbitrary quantity F = F(xi,B) at 
// radial grid point ir and low-field side pitch xi0.
real_t AnalyticBRadialGridGenerator::BounceAverageQty(len_t ir, real_t xi0, "@(xi,B) F(xi,B)"" ){
    real_t theta_b1;
    real_t theta_b2;

    xi = sign(xi0) * sqrt( 1- B/Bmin * (1-xi0^2));

    setBouncePoints(ir,xi0, &theta_b1,&theta_b2); 
    // passing:
    // theta_b1 = -pi
    // theta_b2 = pi
    if (IsTrapped){
        return (1/Vp) * integral( sqrt(g) *(F(xi(theta),B(theta)) + F(-xi(theta),B(theta)) )/2 ,theta, theta_b1, theta_b2  );
    } else {
        return (1/Vp) *  integral( sqrt(g) *(F(xi(theta),B(theta)),theta, theta_b1, theta_b2  );
    }
}
*/

