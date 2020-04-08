/**
 * Implementation of a radial grid in analytic toroidal magnetic geometry 
 * (with given magnetic-axis major radius, plasma minor radius and profiles of elongation, 
 *  triangularity and Shafranov shift as well as reference poloidal flux profile)
 */

#include <cmath>
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <functional>
#include <gsl/gsl_integration.h>

using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * nx: Number of radial grid points.
 * G: Toroidal magnetic field component as function of minor radius 
 * Psi_p0: Reference poloidal magnetic flux as function of minor radius 
 * x0: Value of inner radial flux grid point.
 * xa: Value of outer radial flux grid point.
 */
AnalyticBRadialGridGenerator::AnalyticBRadialGridGenerator(
     const len_t nr,  real_t r0,  real_t ra, real_t R0,
    std::function<real_t(real_t)> G,  std::function<real_t(real_t)> Psi_p0, 
    std::function<real_t(real_t)> kappa, std::function<real_t(real_t)> delta, 
    std::function<real_t(real_t)> Delta
) : RadialGridGenerator(nr), rMin(r0), rMax(ra) {
    this->rMin = r0;
    this->rMax = ra;
    this->R0   = R0;
    this->Btor_G       = G;
    this->Psi_p0       = Psi_p0;
    this->kappa_Elong  = kappa;
    this->delta_triang = delta;
    this->Delta_shafr  = Delta;

}

/*************************************
 * PUBLIC METHODS                    *
 *************************************/
/**
 * (Re-)builds the given radial grid.
 *
 * rGrid: Radial grid to re-build.
 */
bool AnalyticBRadialGridGenerator::Rebuild(const real_t, RadialGrid *rGrid) {
    real_t
        *r    = new real_t[GetNr()],
        *r_f  = new real_t[GetNr()+1],
        *dr   = new real_t[GetNr()],
        *dr_f = new real_t[GetNr()-1];

    // Construct flux grid
    for (len_t i = 0; i < GetNr(); i++)
        dr[i] = (rMax - rMin) / GetNr();

    for (len_t i = 0; i < GetNr()+1; i++)
        r_f[i] = rMin + i*dr[0];

    // Construct cell grid
    for (len_t i = 0; i < GetNr(); i++)
        r[i] = 0.5 * (r_f[i+1] + r_f[i]);

    for (len_t i = 0; i < GetNr()-1; i++)
        dr_f[i] = r[i+1] - r[i];

    rGrid->Initialize(r, r_f, dr, dr_f);


    this->isBuilt = true;

    return true;
}

real_t AnalyticBRadialGridGenerator::diffFunc(real_t r, std::function<real_t(real_t)> F){
    real_t sqrteps = sqrt(__DBL_EPSILON__);
    real_t h = sqrteps * ( 1 + abs(r) ); 
    return (F(r+h)-F(r))/h;
}

/**
 * Re-build the phase space jacobians.
 *
 * grid: Grid to build jacobians for.
 */
void AnalyticBRadialGridGenerator::RebuildJacobians(RadialGrid *rGrid, MomentumGrid **momentumGrids) {
    
    real_t
        **Vp      = new real_t*[GetNr()],
        **Vp_fr   = new real_t*[GetNr()+1],
        **Vp_f1   = new real_t*[GetNr()],
        **Vp_f2   = new real_t*[GetNr()];


    // Construct magnetic field quantities that depend on poloidal angle theta
    real_t kappa, Delta, delta, kappaPrime, DeltaPrime, deltaPrime, psiPrime;
    real_t
        *theta        = new real_t[ntheta],
        *weightsTheta = new real_t[ntheta],
        *st           = new real_t[ntheta],
        *ct           = new real_t[ntheta],
        *R            = new real_t[GetNr()*ntheta],
        *nabla_r2     = new real_t[GetNr()*ntheta],
        *JacobianJ    = new real_t[GetNr()*ntheta],
        *B            = new real_t[GetNr()*ntheta],
        *Bmin         = new real_t[GetNr()],
        *R_f          = new real_t[(GetNr()+1)*ntheta],
        *nabla_r2_f   = new real_t[(GetNr()+1)*ntheta],
        *JacobianJ_f  = new real_t[(GetNr()+1)*ntheta],
        *B_f          = new real_t[(GetNr()+1)*ntheta],
        *Bmin_f       = new real_t[GetNr()+1];

    // Distributing poloidal-angle grid points theta according to the Gauss-Legendre
    // quadrature rule for more efficient flux surface averaging
    const gsl_integration_fixed_type *legendreTypeQuad = gsl_integration_fixed_legendre;
    gsl_integration_fixed_workspace *gsl_w = gsl_integration_fixed_alloc(legendreTypeQuad,ntheta,0,2*M_PI,0,0);
    for(len_t it=0; it<ntheta; it++) {
        theta[it] = gsl_w->x[it];
        weightsTheta[it] = gsl_w->weights[it];
        st[it] = sin(theta[it]);
        ct[it] = cos(theta[it]);
    }
    this->theta = theta;
    this->weightsTheta = weightsTheta;

    const real_t *r;   
    const real_t *r_f; 
    
    for (len_t ir = 0; ir < GetNr(); ir++){
        r   = rGrid->GetR();
        r_f = rGrid->GetR_f();
        kappa = kappa_Elong(r[ir]);
        delta = delta_triang(r[ir]);
        Delta = Delta_shafr(r[ir]);
        kappaPrime = diffFunc(r[ir],kappa_Elong);
        deltaPrime = diffFunc(r[ir],delta_triang);
        DeltaPrime = diffFunc(r[ir],Delta_shafr);
        psiPrime = diffFunc(r[ir],Psi_p0);



        for(len_t it=0; it<ntheta; it++){
            R[ir*ntheta + it] = R0 + Delta + r[ir]*cos(theta[it] + delta*st[it]);
            
            JacobianJ[ir*ntheta+it] = kappa*r[ir]*R[ir*ntheta + it] * ( cos(delta*st[it]) + DeltaPrime*ct[it]
            + st[it]*sin(theta[it]+delta*st[it]) * ( r[ir]*kappaPrime/kappa + delta*ct[it]
            * ( 1 + r[ir]*kappaPrime/kappa - r[ir]*deltaPrime/delta ) ) );
            
            nabla_r2[ir*ntheta + it] = kappa*kappa*r[ir]*r[ir]*R[ir*ntheta+it]*R[ir*ntheta+it]
                *( ct[it]*ct[it] + (1+delta*delta)*(1+delta*delta)/(kappa*kappa) 
                * sin(theta[it]+delta*st[it])*sin(theta[it]+delta*st[it]) ) 
                / ( JacobianJ[ir*ntheta+it]*JacobianJ[ir*ntheta+it] ); 
            
            B[ir*ntheta + it] = Btor_G(r[ir])*Btor_G(r[ir])/(R[ir*ntheta + it]*R[ir*ntheta + it])
                                + nabla_r2[ir*ntheta + it] * psiPrime*psiPrime;
        }    
    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        r   = rGrid->GetR();
        r_f = rGrid->GetR_f();
        kappa = kappa_Elong(r_f[ir]);
        delta = delta_triang(r_f[ir]);
        Delta = Delta_shafr(r_f[ir]);
        kappaPrime = diffFunc(r_f[ir],kappa_Elong);
        deltaPrime = diffFunc(r_f[ir],delta_triang);
        DeltaPrime = diffFunc(r_f[ir],Delta_shafr);
        psiPrime = diffFunc(r_f[ir],Psi_p0);

        for(len_t it=0; it<ntheta; it++){
            R_f[ir*ntheta + it] = R0 + Delta + r_f[ir]*cos(theta[it] + delta*st[it]);
            
            JacobianJ_f[ir*ntheta+it] = kappa*r_f[ir]*R_f[ir*ntheta + it] * ( cos(delta*st[it]) + DeltaPrime*ct[it]
            + st[it]*sin(theta[it]+delta*st[it]) * ( r_f[ir]*kappaPrime/kappa + delta*ct[it]
            * ( 1 + r_f[ir]*kappaPrime/kappa - r_f[ir]*deltaPrime/delta ) ) );
            
            nabla_r2_f[ir*ntheta + it] = kappa*kappa*r_f[ir]*r_f[ir]*R_f[ir*ntheta+it]*R_f[ir*ntheta+it]
                *( ct[it]*ct[it] + (1+delta*delta)*(1+delta*delta)/(kappa*kappa) 
                * sin(theta[it]+delta*st[it])*sin(theta[it]+delta*st[it]) ) 
                / ( JacobianJ[ir*ntheta+it]*JacobianJ[ir*ntheta+it] ); 
            
            B_f[ir*ntheta + it] = Btor_G(r_f[ir])*Btor_G(r_f[ir])/(R_f[ir*ntheta + it]*R_f[ir*ntheta + it])
                                + nabla_r2[ir*ntheta + it] * psiPrime * psiPrime;
        }    
    }
    
    for (len_t ir = 0; ir < GetNr(); ir++){
        Bmin[ir] = B[ir*ntheta];
        for(len_t it=0; it<ntheta; it++)
            if (Bmin[ir] > B[ir*ntheta+it])
                Bmin[ir] = B[ir*ntheta+it];
        }

    for (len_t ir = 0; ir < GetNr()+1; ir++){
        Bmin_f[ir] = B_f[ir*ntheta];
        for(len_t it=0; it<ntheta; it++)
            if (Bmin_f[ir] > B_f[ir*ntheta+it])
                Bmin_f[ir] = B_f[ir*ntheta+it];
        }

    
    
    rGrid->InitializeMagneticField(
        ntheta, theta, B, B_f, Bmin, Bmin_f, JacobianJ, JacobianJ_f
    );

    
    // Set Vp, Vp_f1 and Vp_f2
    for (len_t ir = 0; ir < GetNr(); ir++) {
        const MomentumGrid *mg = momentumGrids[ir];
        //const real_t r = rGrid->GetR(ir);

        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();

        Vp[ir]    = new real_t[n1*n2];
        Vp_f1[ir] = new real_t[(n1+1)*n2];
        Vp_f2[ir] = new real_t[n1*(n2+1)];
        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                Vp[ir][j*n1 + i] = EvaluateBounceSurfaceIntegral(rGrid, mg, ir, i, j, 0,[&](real_t,real_t){return 1;});;
            }
        }

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1+1; i++) {
                Vp_f1[ir][j*(n1+1) + i] = EvaluateBounceSurfaceIntegral(rGrid, mg, ir, i, j, 2,[&](real_t,real_t){return 1;});;
            }
        }

        for (len_t j = 0; j < n2+1; j++) {
            for (len_t i = 0; i < n1; i++) {
                Vp_f2[ir][j*n1 + i] = EvaluateBounceSurfaceIntegral(rGrid, mg, ir, i, j, 3,[&](real_t,real_t){return 1;});
            }
        }
    }


    // Set Vp_fr
    for (len_t ir = 0; ir < GetNr()+1; ir++) {
        // XXX: We inherently assume that the momentum grids at all
        // radii are the same here, so we might as well just evaluate
        // everything on the innermost momentum grid
        const MomentumGrid *mg = momentumGrids[0];

        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();

        Vp_fr[ir] = new real_t[n1*n2];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                Vp_fr[ir][j*n1 + i] = EvaluateBounceSurfaceIntegral(rGrid, mg, ir, i, j, 1,[&](real_t,real_t){return 1;});
            }
        }
    }

    rGrid->InitializeVprime(Vp, Vp_fr, Vp_f1, Vp_f2);
    delete [] st;
    delete [] ct;
}

/**
 * Re-build flux surface averaged quantities.
 *
 * grid: Grid to build them on.
 */
void AnalyticBRadialGridGenerator::RebuildFSAvgQuantities(RadialGrid *rGrid, MomentumGrid **momentumGrids) {
    real_t
     *effectivePassingFraction = new real_t[GetNr()],
     *magneticFieldMRS         = new real_t[GetNr()],
     *nabla_r2_avg             = new real_t[GetNr()],
     **xiBounceAverage_f1      = new real_t*[GetNr()],
     **xiBounceAverage_f2      = new real_t*[GetNr()],
     **xi21MinusXi2OverB2_f1   = new real_t*[GetNr()],
     **xi21MinusXi2OverB2_f2   = new real_t*[GetNr()];
     

    for (len_t ir = 0; ir < GetNr(); ir++) {
        // TODO
        effectivePassingFraction[ir] = 1; // = 4/3 * integral( x/ FluxSurfaceAverage(ir, sqrt(1-x*B/Bmax),x,0,1 ) )
        magneticFieldMRS[ir] = 1; // = sqrt( FluxSurfaceAverage(ir, B*B )) 

        const MomentumGrid *mg = momentumGrids[ir];
        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();
        
        xiBounceAverage_f1[ir] = new real_t[(n1+1)*n2];
        xiBounceAverage_f2[ir] = new real_t[n1*(n2+1)];
        xi21MinusXi2OverB2_f1[ir] = new real_t[(n1+1)*n2];
        xi21MinusXi2OverB2_f2[ir] = new real_t[n1*(n2+1)];
         
//        const real_t *p1     = mg->GetP1();
//        const real_t *p1_f   = mg->GetP1_f();
//        const real_t *p2     = mg->GetP2();
//        const real_t *p2_f   = mg->GetP2_f();
        
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

    rGrid->InitializeFSAvg(effectivePassingFraction, magneticFieldMRS,
                            xiBounceAverage_f1, xiBounceAverage_f2,
                            xi21MinusXi2OverB2_f1, xi21MinusXi2OverB2_f2,
                            nabla_r2_avg);
    }
}





/**
 * Calculates the bounce average {F} of an arbitrary function F=F(B/Bmin).
 */
real_t AnalyticBRadialGridGenerator::FluxSurfaceAverageQuantity(RadialGrid *rGrid, len_t ir, bool rFluxGrid, std::function<real_t(real_t)> F){
    

    const real_t *JacobianJ = (
        rFluxGrid ?
            rGrid->GetJacobian_f(ir) :
            rGrid->GetJacobian(ir) 
    );
    const real_t *B = (
        rFluxGrid ?
            rGrid->BOfTheta_f(ir) :
            rGrid->BOfTheta(ir) 
    );
    const real_t Bmin = (
        rFluxGrid ? 
            rGrid->GetBmin_f(ir) :
            rGrid->GetBmin(ir)
    );
    const real_t VolVp = (
        rFluxGrid ?
            rGrid->GetVolVp_f(ir) :
            rGrid->GetVolVp(ir)
    );

    real_t dTh = theta[1]-theta[0];
//    real_t *Vp = rGrid->GetVp(ir);
    real_t FSA;
    for (len_t it=0; it<rGrid->GetNTheta(); it++) {
        FSA += dTh * JacobianJ[it] * F( B[it]/Bmin );
    }

    return FSA / VolVp;
}




real_t AnalyticBRadialGridGenerator::EvaluateBounceSurfaceIntegral(RadialGrid *rGrid, const MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGrid, std::function<real_t(real_t,real_t)> F){

    real_t xi0 = mg->GetXi0(i,j);
    const real_t *B = rGrid->BOfTheta(ir);
    real_t Bmin = rGrid->GetBmin(ir);
    const real_t *Jacobian = rGrid->GetJacobian(ir); 
    real_t p1 = mg->GetP1(i);
    real_t p2 = mg->GetP2(j);

    // if evaluating bounce average on flux grid, update values of certain variables to their flux-grid values
    switch (fluxGrid){ 
        case 1: // r flux grid
            B = rGrid->BOfTheta_f();
            Bmin = rGrid->GetBmin_f(ir);
            Jacobian = rGrid->GetJacobian_f(ir);
            break;
        case 2: // p1 flux grid
            xi0 = mg->GetXi0_f1(i,j);
            p1 = mg->GetP1_f(i);
            break;
        case 3: // p2 flux grid
            xi0 = mg->GetXi0_f2(i,j);
            p2 = mg->GetP2_f(j);
            break;
      
    }


    real_t *sqrtg = new real_t[ntheta];    
    mg->EvaluateMetric(p1, p2, ir, rGrid, ntheta, theta, fluxGrid==1, sqrtg);

    /*
    // The following block finds the indices it_bounce[0] and it_bounce[1] corresponding to 
    // the two bounce points theta_bounce1 and theta_bounce2
    bool sign_prev = (1-Bmin/B[0]*(1-xi0*xi0) >= 0); // true if sign is positive
    bool sign_this;
    len_t offset;
    len_t it_bounce[2]; 
    len_t store_ind=0;
    real_t xi2;
    for (len_t it=1; it<ntheta; it++){
        xi2 = 1-Bmin/B[it]*(1-xi0*xi0);     // particle does not reach theta[it] if xi2<0 
        sign_this = (xi2 >=0 );
        if (sign_this != sign_prev){        // if sign changes, store index to it_bounce
            if (sign_this)
                it_bounce[1] = it;
            else
                it_bounce[0] = it-1;       // if we are now in the trapped region, store previous index
            it_bounce[store_ind] = it-offset;
            store_ind += 1;
        }
    } 
    */
    bool *onOrbit = new bool[ntheta]; // true for all theta between theta_bounce1 and theta_bounce2
    bool isTrapped = false;
    for (len_t it=0; it<ntheta; it++){
        onOrbit[it] = (1-Bmin/B[it]*(1-xi0*xi0) > 0);
        if (!onOrbit[it])
            isTrapped = true; //if there is any 
    }
    
    std::function<real_t(real_t,real_t)> F_eff;
    
    if (isTrapped)
        F_eff = [&](real_t x, real_t  y){return (F(x,y) + F(-x,y))/2;};
    else 
        F_eff = F;

    real_t sign;
    if (xi0>=0) 
        sign=1; 
    else 
        sign=-1;

    real_t xi;
    real_t bounceIntegral = 0;
    real_t *trapzWeights = new real_t[ntheta];
    trapzWeights[0] = (theta[1]-theta[0])/2;
    trapzWeights[ntheta-1] = (theta[ntheta-1]-theta[ntheta-2])/2;
    for (len_t it=1; it<ntheta-1; it++)
        trapzWeights[it] = theta[it] - theta[it-1];

    // trapezoidal quadrature for passing orbits, riemann for trapped. 
    // Should interpolate to new theta grid for accuracy.
    for (len_t it=0; it<ntheta; it++){
        if (onOrbit[it]){
            xi = sign*sqrt(1-Bmin/B[it]*(1-xi0*xi0));
            bounceIntegral += trapzWeights[it]*Jacobian[it]*sqrtg[it]*F_eff(xi,B[it]/Bmin);
        }

    }

    return bounceIntegral;

}
/**
 * Calculates the bounce average {F} of an arbitrary function F=F(xi,B/Bmin).
 * We may want to improve this by interpolating quantities to a finer where, perhaps, 
 * there are grid points on +/- theta_bounce
 */
real_t AnalyticBRadialGridGenerator::BounceAverageQuantity(RadialGrid *rGrid, const MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGrid, std::function<real_t(real_t,real_t)> F){
    // FluxGrid: 0: distribution grid. 1: r flux grid. 2: p1 flux grid. 3: p2 flux grid.

    len_t np1 = mg->GetNp1();
    real_t Vp;
    switch (fluxGrid){
        case 0:
            Vp = rGrid->GetVp(ir)[j*np1+i];
            break;
        case 1:
            Vp = rGrid->GetVp_fr(ir)[j*np1+i];
            break;
        case 2:
            Vp = rGrid->GetVp_f1(ir)[j*(np1+1)+i];
            break;
        case 3:
            Vp = rGrid->GetVp_f2(ir)[j*np1+1];
            break;
        default:
            Vp = 0; // return error?
    }

    
    return EvaluateBounceSurfaceIntegral(rGrid,mg,ir,i,j,fluxGrid, F) / Vp;

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

