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
#include <gsl/gsl_spline.h>

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
    real_t *rProfiles, len_t nrProfiles, real_t *Gs, real_t *psi_p0s,
             real_t *kappas, real_t *deltas, real_t *Deltas
) : RadialGridGenerator(nr), rMin(r0), rMax(ra), R0(R0){
    
    this->GsProvided     = Gs;
    this->psisProvided   = psi_p0s;
    this->kappasProvided = kappas;
    this->deltasProvided = deltas;
    this->DeltasProvided = Deltas;
    this->nrProfiles     = nrProfiles;
    this->rProfilesProvided = rProfiles;
    
    spline_x = gsl_spline_alloc(gsl_interp_steffen, nrProfiles);
    gsl_acc  = gsl_interp_accel_alloc(); 
}

AnalyticBRadialGridGenerator::~AnalyticBRadialGridGenerator(){
    gsl_spline_free (spline_x);
    gsl_interp_accel_free (gsl_acc);
}


/**
 * Interpolates input shape-parameter profiles (kappa, delta, ...) which are defined on 
 * input rProfilesProvided array to the r and r_f grids
 */
void AnalyticBRadialGridGenerator::InterpolateInputProfileToGrid(real_t *r, real_t *r_f, real_t *x,real_t *xPrime, real_t *x_f, real_t *xPrime_f,real_t *xProvided){
    x        = new real_t[GetNr()];
    xPrime   = new real_t[GetNr()];
    x_f      = new real_t[GetNr()+1];
    xPrime_f = new real_t[GetNr()+1];

    gsl_spline_init(spline_x, rProfilesProvided, xProvided, nrProfiles);
    real_t sqrteps = sqrt(__DBL_EPSILON__);
    real_t h;
    for (len_t ir=0; ir<GetNr(); ir++){
        x[ir]      = gsl_spline_eval(spline_x,r[ir],gsl_acc);
        
        h = sqrteps * ( 1 + abs(r[ir]) );
        if (ir==GetNr()-1)
            h = -h;
        xPrime[ir] = (gsl_spline_eval(spline_x,r[ir]+h,gsl_acc) - x[ir])/h;
        }
    for (len_t ir=0; ir<GetNr()+1; ir++){
        x_f[ir]      = gsl_spline_eval(spline_x,r_f[ir],gsl_acc);
        
        h = sqrteps * ( 1 + abs(r_f[ir]) );
        if (ir==GetNr())
            h = -h;
        xPrime_f[ir] = (gsl_spline_eval(spline_x,r[ir]+h,gsl_acc) - x_f[ir])/h;
        }
}


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

    InterpolateInputProfileToGrid(r,r_f,G,GPrime,G_f,GPrime_f,GsProvided);
    InterpolateInputProfileToGrid(r,r_f,psi,psiPrime,psi_f,psiPrime_f,psisProvided);
    InterpolateInputProfileToGrid(r,r_f,kappa,kappaPrime,kappa_f,kappaPrime_f,kappasProvided);
    InterpolateInputProfileToGrid(r,r_f,delta,deltaPrime,delta_f,deltaPrime_f,deltasProvided);
    InterpolateInputProfileToGrid(r,r_f,Delta,DeltaPrime,Delta_f,DeltaPrime_f,DeltasProvided);
    
    return true;
}

real_t AnalyticBRadialGridGenerator::diffFunc(real_t r, std::function<real_t(real_t)> F){
    real_t sqrteps = sqrt(__DBL_EPSILON__);
    real_t h = sqrteps * ( 1 + abs(r) ); 
    return (F(r+h/2)-F(r-h/2))/h;
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
    // real_t kappa, Delta, delta, kappaPrime, DeltaPrime, deltaPrime, psiPrime;
    real_t
        *theta        = new real_t[ntheta],
        *weightsTheta = new real_t[ntheta],
        *st           = new real_t[ntheta],
        *ct           = new real_t[ntheta],
        *JacobianJ    = new real_t[GetNr()*ntheta],
        *B            = new real_t[GetNr()*ntheta],
        *Bmin         = new real_t[GetNr()],
        *JacobianJ_f  = new real_t[(GetNr()+1)*ntheta],
        *B_f          = new real_t[(GetNr()+1)*ntheta],
        *Bmin_f       = new real_t[GetNr()+1];

    R          = new real_t[GetNr()*ntheta];
    R_f        = new real_t[(GetNr()+1)*ntheta];
    nabla_r2   = new real_t[GetNr()*ntheta];
    nabla_r2_f = new real_t[(GetNr()+1)*ntheta];
        
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

        for(len_t it=0; it<ntheta; it++){
            R[ir*ntheta + it] = R0 + Delta[ir] + r[ir]*cos(theta[it] + delta[ir]*st[it]);
            
            JacobianJ[ir*ntheta+it] = kappa[ir]*r[ir]*R[ir*ntheta + it] * ( cos(delta[ir]*st[it]) + DeltaPrime[ir]*ct[it]
            + st[it]*sin(theta[it]+delta[ir]*st[it]) * ( r[ir]*kappaPrime[ir]/kappa[ir] + delta[ir]*ct[it]
            * ( 1 + r[ir]*kappaPrime[ir]/kappa[ir] - r[ir]*deltaPrime[ir]/delta[ir] ) ) );
            
            nabla_r2[ir*ntheta + it] = kappa[ir]*kappa[ir]*r[ir]*r[ir]*R[ir*ntheta+it]*R[ir*ntheta+it]
                *( ct[it]*ct[it] + (1+delta[ir]*delta[ir])*(1+delta[ir]*delta[ir])/(kappa[ir]*kappa[ir]) 
                * sin(theta[it]+delta[ir]*st[it])*sin(theta[it]+delta[ir]*st[it]) ) 
                / ( JacobianJ[ir*ntheta+it]*JacobianJ[ir*ntheta+it] ); 
            
            B[ir*ntheta + it] = G[ir]*G[ir]/(R[ir*ntheta + it]*R[ir*ntheta + it])
                                + nabla_r2[ir*ntheta + it] * psiPrime[ir]*psiPrime[ir];
        }    
    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        r   = rGrid->GetR();
        r_f = rGrid->GetR_f();

        for(len_t it=0; it<ntheta; it++){
            R_f[ir*ntheta + it] = R0 + Delta_f[ir] + r_f[ir]*cos(theta[it] + delta_f[ir]*st[it]);
            
            JacobianJ_f[ir*ntheta+it] = kappa_f[ir]*r_f[ir]*R_f[ir*ntheta + it] * ( cos(delta_f[ir]*st[it]) + DeltaPrime_f[ir]*ct[it]
            + st[it]*sin(theta[it]+delta_f[ir]*st[it]) * ( r_f[ir]*kappaPrime_f[ir]/kappa_f[ir] + delta_f[ir]*ct[it]
            * ( 1 + r_f[ir]*kappaPrime_f[ir]/kappa_f[ir] - r_f[ir]*deltaPrime_f[ir]/delta_f[ir] ) ) );
            
            nabla_r2_f[ir*ntheta + it] = kappa_f[ir]*kappa_f[ir]*r_f[ir]*r_f[ir]*R_f[ir*ntheta+it]*R_f[ir*ntheta+it]
                *( ct[it]*ct[it] + (1+delta_f[ir]*delta_f[ir])*(1+delta_f[ir]*delta_f[ir])/(kappa_f[ir]*kappa_f[ir]) 
                * sin(theta[it]+delta_f[ir]*st[it])*sin(theta[it]+delta_f[ir]*st[it]) ) 
                / ( JacobianJ[ir*ntheta+it]*JacobianJ[ir*ntheta+it] ); 
            
            B_f[ir*ntheta + it] = G_f[ir]*G_f[ir]/(R_f[ir*ntheta + it]*R_f[ir*ntheta + it])
                                + nabla_r2[ir*ntheta + it] * psiPrime_f[ir] * psiPrime_f[ir];
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
        ntheta, theta, B, B_f, Bmin, Bmin_f, JacobianJ, JacobianJ_f);


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
 * The function is used in the evaluation of the effective passing fraction,
 * and represents x / <1-x B/Bmax>
 */
struct EPF_params {real_t BminOverBmax; len_t ir;  AnalyticBRadialGridGenerator *ABrgg; RadialGrid *rg; };
real_t AnalyticBRadialGridGenerator::effectivePassingFractionIntegrand(real_t x, void *p){
    struct EPF_params *params = (struct EPF_params *) p;
    AnalyticBRadialGridGenerator *ABrgg = params->ABrgg;
    RadialGrid *rg = params->rg;
    real_t BminOverBmax = params->BminOverBmax; 
    len_t ir = params->ir;
    std::function<real_t(real_t)> fluxAvgFunc = [x,BminOverBmax](real_t BOverBmin){
        return sqrt(1 - x * BminOverBmax * BOverBmin );
    };
    return x/ ABrgg->FluxSurfaceAverageQuantity(rg, ir, false, fluxAvgFunc);
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
     *magneticFieldMRS_f       = new real_t[GetNr()+1],
     *nablaR2OverR2_avg        = new real_t[GetNr()],
     *nablaR2OverR2_avg_f      = new real_t[GetNr()+1],
     *OneOverR2_avg            = new real_t[GetNr()],
     *OneOverR2_avg_f          = new real_t[GetNr()+1],
     **xiBounceAverage_f1      = new real_t*[GetNr()],
     **xiBounceAverage_f2      = new real_t*[GetNr()],
     **xi21MinusXi2OverB2_f1   = new real_t*[GetNr()],
     **xi21MinusXi2OverB2_f2   = new real_t*[GetNr()],
     **OneOverBOverXi_avg_f1   = new real_t*[GetNr()],
     **OneOverBOverXi_avg_f2   = new real_t*[GetNr()];
        

    real_t *nablaR2OverR2 = new real_t[ntheta]; 
    real_t *OneOverR2     = new real_t[ntheta]; 

    real_t Bmin, Bmax, EPF_integral;
    real_t xi0_f1,xi0_f2,BOverXi_avg,sign_xi0;
    bool isTrapped;
    gsl_integration_workspace *gsl_w = gsl_integration_workspace_alloc(1000);
    for (len_t ir = 0; ir < GetNr(); ir++) {
        magneticFieldMRS[ir] = rGrid->GetBmin(ir)*sqrt( FluxSurfaceAverageQuantity(rGrid,ir,false,[](real_t BOverBmin){return BOverBmin*BOverBmin;}) ) ;

        for (len_t it=0; it<ntheta; it++){
            nablaR2OverR2[it] = R0*R0*nabla_r2[ir*ntheta+it]/(R[ir*ntheta+it]*R[ir*ntheta+it]);
            OneOverR2[it]     = R0*R0/(R[ir*ntheta+it]*R[ir*ntheta+it]);
        }

        nablaR2OverR2_avg[ir] = FluxSurfaceAverageQuantity(rGrid,ir,false,nablaR2OverR2 );
        OneOverR2_avg[ir] = FluxSurfaceAverageQuantity(rGrid,ir,false,OneOverR2 );

        /**
         * The following block integrates the effective passing fraction.
         */
        gsl_function EPF_func;
        const real_t *B = rGrid->BOfTheta(ir);
        Bmax = B[0];
        for (len_t it = 0; it<ntheta; it++){
            if(Bmax<B[it])
                Bmax = B[it];
        }
        Bmin = rGrid->GetBmin(ir);
        real_t BminOverBmax = Bmin/Bmax;
        EPF_params paramstruct = {BminOverBmax,ir,this,rGrid}; 
        EPF_func.function = &(effectivePassingFractionIntegrand);
        EPF_func.params = &paramstruct;
        gsl_integration_qags(&EPF_func, 0,1,0,1e-7,1000,gsl_w,&EPF_integral,nullptr );
        effectivePassingFraction[ir] = (3/4) * magneticFieldMRS[ir]*magneticFieldMRS[ir]  / (Bmax*Bmax) * EPF_integral;
        

        const MomentumGrid *mg = momentumGrids[ir];
        const len_t n1 = mg->GetNp1();
        const len_t n2 = mg->GetNp2();
        
        xiBounceAverage_f1[ir]     = new real_t[(n1+1)*n2];
        xiBounceAverage_f2[ir]     = new real_t[n1*(n2+1)];
        xi21MinusXi2OverB2_f1[ir]  = new real_t[(n1+1)*n2];
        xi21MinusXi2OverB2_f2[ir]  = new real_t[n1*(n2+1)];
        OneOverBOverXi_avg_f1[ir]  = new real_t[(n1+1)*n2];
        OneOverBOverXi_avg_f2[ir]  = new real_t[n1*(n2+1)];
        
        
        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1+1; i++) {
                xiBounceAverage_f1[ir][j*(n1+1)+i]    = BounceAverageQuantity(rGrid, mg, ir, i, j, 2, [](real_t xi, real_t  ){return xi;} );
                xi21MinusXi2OverB2_f1[ir][j*(n1+1)+i] = BounceAverageQuantity(rGrid, mg, ir, i, j, 2, [](real_t xi, real_t BOverBMin ){return xi*xi*(1-xi*xi)/(BOverBMin*BOverBMin);} );
                
                xi0_f1 = mg->GetXi0_f1(i,j);
                sign_xi0 = (xi0_f1 > 0) - (xi0_f1 <0) ;
                isTrapped =  (1-xi0_f1*xi0_f1)/BminOverBmax > 1 ;
                OneOverBOverXi_avg_f1[ir][j*(n1+1)+i] = 0;
                if (!isTrapped){
                    if (!(Bmax==Bmin && xi0_f1==0)){
                        BOverXi_avg = Bmin * FluxSurfaceAverageQuantity(rGrid, ir, false, [xi0_f1](real_t BoBm){return BoBm/sqrt(1-BoBm*(1-xi0_f1*xi0_f1));});
                        OneOverBOverXi_avg_f1[ir][j*(n1+1)+i] +=  sign_xi0 * magneticFieldMRS[ir]/BOverXi_avg;
                    }
                }
            }
        }
        for (len_t j = 0; j < n2+1; j++) {
            for (len_t i = 0; i < n1; i++) {
                xiBounceAverage_f2[ir][j*n1+i]    = BounceAverageQuantity(rGrid, mg, ir, i, j, 3, [](real_t xi, real_t  ){return xi;} );
                xi21MinusXi2OverB2_f2[ir][j*n1+i] = BounceAverageQuantity(rGrid, mg, ir, i, j, 3, [](real_t xi, real_t BOverBMin ){return xi*xi*(1-xi*xi)/(BOverBMin*BOverBMin);} );

                xi0_f2 = mg->GetXi0_f2(i,j);
                sign_xi0 = (xi0_f2 > 0) - (xi0_f2 <0) ;
                isTrapped =  (1-xi0_f2*xi0_f2)/BminOverBmax > 1 ;
                OneOverBOverXi_avg_f2[ir][j*(n1+1)+i] = 0;
                if (!isTrapped){
                    if (!(Bmax==Bmin && xi0_f2==0)){
                        BOverXi_avg = Bmin * FluxSurfaceAverageQuantity(rGrid, ir, false, [xi0_f2](real_t BoBm){return BoBm/sqrt(1-BoBm*(1-xi0_f2*xi0_f2));});
                        OneOverBOverXi_avg_f2[ir][j*n1+i] +=  sign_xi0 * magneticFieldMRS[ir]/BOverXi_avg;
                    }
                }
            }
        }
    }

    real_t *nablaR2OverR2_f = new real_t[ntheta]; 
    real_t *OneOverR2_f = new real_t[ntheta]; 
    for (len_t ir = 0; ir < GetNr()+1; ir++) {
        magneticFieldMRS_f[ir] = rGrid->GetBmin_f(ir)*sqrt( FluxSurfaceAverageQuantity(rGrid,ir,true,[](real_t BOverBmin){return BOverBmin*BOverBmin;}) ) ;
        for (len_t it=0; it<ntheta; it++) {
            nablaR2OverR2_f[it] = R0*R0*nabla_r2_f[ir*ntheta+it]/(R_f[ir*ntheta+it]*R_f[ir*ntheta+it]);
            OneOverR2_f[it]     = R0*R0/(R_f[ir*ntheta+it]*R_f[ir*ntheta+it]);
        }
        nablaR2OverR2_avg_f[ir] = FluxSurfaceAverageQuantity(rGrid,ir,true,nablaR2OverR2_f );
        OneOverR2_avg_f[ir] = FluxSurfaceAverageQuantity(rGrid,ir,true,OneOverR2_f );
    }
    rGrid->InitializeFSAvg(effectivePassingFraction, magneticFieldMRS,magneticFieldMRS_f,
                            xiBounceAverage_f1, xiBounceAverage_f2,
                            xi21MinusXi2OverB2_f1, xi21MinusXi2OverB2_f2,
                            nablaR2OverR2_avg, nablaR2OverR2_avg_f,
                            OneOverR2_avg, OneOverR2_avg_f,
                            OneOverBOverXi_avg_f1,OneOverBOverXi_avg_f2);

    gsl_integration_workspace_free (gsl_w);
    delete [] nablaR2OverR2;
    delete [] nablaR2OverR2_f;
    delete [] nabla_r2;
    delete [] nabla_r2_f;
    delete [] R;
    delete [] R_f;
}





/**
 * Calculates the flux surface integral V'<F> of an arbitrary function F=F(B/Bmin).
 */
real_t AnalyticBRadialGridGenerator::EvaluateFluxSurfaceIntegral(RadialGrid *rGrid, len_t ir, bool rFluxGrid, std::function<real_t(real_t)> F){
    const real_t *B;
    real_t Bmin;
    if(rFluxGrid){
        B = rGrid->BOfTheta_f(ir);
        Bmin = rGrid->GetBmin_f(ir);
    } else {
        B = rGrid->BOfTheta(ir);
        Bmin = rGrid->GetBmin(ir);
    }

    // Does this initialization cause memory issues? where do all the F_arrays go..?
    real_t *F_array = new real_t[ntheta]; 
    for (len_t it=0; it<ntheta; it++) {
        F_array[it] = F( B[it]/Bmin );
    }

    return EvaluateFluxSurfaceIntegral(rGrid,ir,rFluxGrid,F);
}


/**
 * Calculates the flux surface integral of an array F on theta. (size ntheta)
 */
real_t AnalyticBRadialGridGenerator::EvaluateFluxSurfaceIntegral(RadialGrid *rGrid, len_t ir, bool rFluxGrid, real_t *F){
    const real_t *Jacobian;
    if(rFluxGrid){
        Jacobian = rGrid->GetJacobian_f(ir);
    } else {
        Jacobian = rGrid->GetJacobian(ir);
    }

    real_t fluxIntegral = 0;
    for (len_t it=0; it<rGrid->GetNTheta(); it++) {
        fluxIntegral += weightsTheta[it] * Jacobian[it] * F[it];
    }
    return fluxIntegral;
}




/**
 * Integrates an arbitrary function F(xi,B/Bmin) over sqrt(g) dtheta
 * between the bounce points, with sqrt(g) the phase-space jacobian
 */
real_t AnalyticBRadialGridGenerator::EvaluateBounceSurfaceIntegral(RadialGrid *rGrid, const MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGrid, std::function<real_t(real_t,real_t)> F){
    real_t *sqrtg   = new real_t[ntheta];    
    bool *onOrbit   = new bool[ntheta];   // true for all theta between theta_bounce1 and theta_bounce2
    real_t *weights = new real_t[ntheta]; // quadrature for integral
    
    
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

    mg->EvaluateMetric(p1, p2, ir, rGrid, ntheta, theta, fluxGrid==1, sqrtg);

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

    

    if (isTrapped){
        weights[0] = (theta[1]-theta[0])/2;
        weights[ntheta-1] = (theta[ntheta-1]-theta[ntheta-2])/2;
        for (len_t it=1; it<ntheta-1; it++)
            weights[it] = theta[it] - theta[it-1];
    } else {
        for (len_t it=0; it<ntheta; it++)
            weights[it] = weightsTheta[it];
    }
    
    real_t sign;
    if (xi0>=0) 
        sign=1; 
    else 
        sign=-1;
    real_t xi;
    real_t bounceIntegral = 0;
    // cool Gauss-Legendre quadrature for passing orbits, crappy riemann sum for trapped. 
    // Should interpolate to new theta grid including bounce points for accuracy.
    for (len_t it=0; it<ntheta; it++){
        if (onOrbit[it]){
            xi = sign*sqrt(1-Bmin/B[it]*(1-xi0*xi0));
            bounceIntegral += weights[it]*Jacobian[it]*sqrtg[it]*F_eff(xi,B[it]/Bmin);
        }

    }
    delete [] weights;
    delete [] onOrbit;
    delete [] sqrtg;
    return bounceIntegral;

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

