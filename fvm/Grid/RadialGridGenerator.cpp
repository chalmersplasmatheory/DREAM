
#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>


using namespace std;
using namespace DREAM::FVM;

RadialGridGenerator::RadialGridGenerator(const len_t nr) : nr(nr) {

}


RadialGridGenerator::~RadialGridGenerator(){
    // DeallocateGridQuantities(momentumGrids); err.. not sure how to do this properly.
    DeallocateMagneticQuantities();
    DeallocateInterpolators();
    DeallocateMagneticFieldData();

}

void RadialGridGenerator::RebuildJacobians(RadialGrid *rGrid, MomentumGrid **momentumGrids) {
    // move this and rGrid->Initialize.. to constructor?
    DeallocateMagneticFieldData();
    CreateMagneticFieldData(rGrid->GetR(),rGrid->GetR_f());

    rGrid->InitializeMagneticField(
        ntheta_ref, theta_ref, B_ref, B_ref_f, Bmin, Bmin_f, Bmax, Bmax_f, Gtor, Gtor_f
    );
    
    InitializeBounceAverage(momentumGrids);

    rGrid->InitializeVprime(GetVp(0),GetVp(1),
                            GetVp(2),GetVp(3),
                            GetVpVol(false), GetVpVol(true));


    DeallocateMagneticQuantities();
}


void RadialGridGenerator::InitializeBounceAverage(MomentumGrid **momentumGrids){

    InitializeMagneticQuantities();

    // if ntheta_ref = 1, we assume cylindrical geometry and will set
    // bounce and flux surface averaging to be identity functions
    // and avoid defining a bunch of stuff.
    if(ntheta_ref==1){
        ntheta_interp=1;
        theta[0]   = 0;
        weights[0] = 1;
        for (len_t ir=0;ir<nr;ir++){
            Bmin[ir] = B_ref[ir][0];
            Bmax[ir] = B_ref[ir][0];
            B[ir][0] = B_ref[ir][0];
            ROverR0[ir][0]  = ROverR0_ref[ir][0];
            Jacobian[ir][0] = Jacobian_ref[ir][0];
            NablaR2[ir][0]  = NablaR2_ref[ir][0];
            
        }
        for (len_t ir=0;ir<nr+1;ir++){
            Bmin_f[ir] = B_ref_f[ir][0];
            Bmax_f[ir] = B_ref_f[ir][0];
            B_f[ir][0] = B_ref_f[ir][0];
            ROverR0_f[ir][0]  = ROverR0_ref_f[ir][0];
            Jacobian_f[ir][0] = Jacobian_ref_f[ir][0];
            NablaR2_f[ir][0]  = NablaR2_ref_f[ir][0];
        }

    } else {

        InitializeInterpolators();

        // Create reference Gauss-Legendre quadrature on x\in[0,1]
        gsl_integration_fixed_workspace *gsl_GL = gsl_integration_fixed_alloc(thetaGridType,ntheta_interp,0,1,0,0);
        this->x_GL_ref = gsl_GL->x;
        this->weights_GL_ref = gsl_GL->weights;

        for (len_t it=0; it<ntheta_interp; it++) {
            theta[it]   = 2*M_PI * x_GL_ref[it];
            weights[it] = 2*M_PI * weights_GL_ref[it];
        }

        
        for (len_t ir=0; ir<nr;ir++){
            Bmin[ir] = B_ref[ir][0];
            Bmax[ir] = B_ref[ir][0];
            for (len_t it = 1; it<ntheta_ref; it++){
                if (Bmin[ir] >  B_ref[ir][it])
                    Bmin[ir] =  B_ref[ir][it];
                if (Bmax[ir] <= B_ref[ir][it])
                    Bmax[ir] =  B_ref[ir][it];
            }
            for (len_t it=0; it<ntheta_interp; it++){
                B[ir][it]        = gsl_spline_eval(B_interpolator[ir], theta[it], gsl_acc);
                ROverR0[ir][it]  = gsl_spline_eval(ROverR0_interpolator[ir], theta[it], gsl_acc);
                Jacobian[ir][it] = gsl_spline_eval(Jacobian_interpolator[ir], theta[it], gsl_acc);
                NablaR2[ir][it]  = gsl_spline_eval(NablaR2_interpolator[ir], theta[it], gsl_acc);
            }
        }

        for (len_t ir=0; ir<nr+1;ir++){
            Bmin_f[ir] = B_ref_f[ir][0];
            Bmax_f[ir] = B_ref_f[ir][0];
            for (len_t it = 1; it<ntheta_ref; it++){
                if (Bmin_f[ir] > B_ref_f[ir][it])
                    Bmin_f[ir] = B_ref_f[ir][it];
                if (Bmax_f[ir] <= B_ref_f[ir][it])
                    Bmax_f[ir] = B_ref_f[ir][it];
            }
            for (len_t it=0; it<ntheta_interp; it++){
                B_f[ir][it]        = gsl_spline_eval(B_interpolator_fr[ir], theta[it], gsl_acc);
                ROverR0_f[ir][it]  = gsl_spline_eval(ROverR0_interpolator_fr[ir], theta[it], gsl_acc);
                NablaR2_f[ir][it]  = gsl_spline_eval(NablaR2_interpolator_fr[ir], theta[it], gsl_acc);
                Jacobian_f[ir][it] = gsl_spline_eval(Jacobian_interpolator_fr[ir], theta[it], gsl_acc);
            }
        }

        EvaluateGrids(momentumGrids);
        
        DeallocateInterpolators();

    }
}





// Evaluates the flux surface average <F> of a function F = F(B/Bmin, R/R0, |nabla r|^2) on radial grid point ir. 
real_t RadialGridGenerator::CalculateFluxSurfaceAverage(len_t ir, bool rFluxGrid, std::function<real_t(real_t,real_t,real_t)> F){
    if (ntheta_interp == 1){
        return F(1,1,1);
    } else 
        return EvaluateFluxSurfaceIntegral(ir,rFluxGrid, F) / GetVpVol(ir,rFluxGrid);
}


real_t RadialGridGenerator::EvaluateFluxSurfaceIntegral(len_t ir, bool rFluxGrid, std::function<real_t(real_t,real_t,real_t)> F){
    real_t Bmin;
    real_t *B, *Jacobian, *NablaR2, *ROverR0;
    
    if (rFluxGrid){
        Bmin     = this->Bmin_f[ir];
        B        = this->B_f[ir]; 
        Jacobian = this->Jacobian_f[ir];
        NablaR2  = this->NablaR2_f[ir];
        ROverR0  = this->ROverR0_f[ir];
    } else {
        Bmin     = this->Bmin[ir];
        B        = this->B[ir];
        Jacobian = this->Jacobian[ir];
        NablaR2  = this->NablaR2[ir];
        ROverR0  = this->ROverR0[ir];
    }

    real_t fluxSurfaceIntegral = 0;

    for (len_t it = 0; it<ntheta_interp; it++){
        fluxSurfaceIntegral += weights[it] * Jacobian[it] * F(B[it]/Bmin, ROverR0[it], NablaR2[it]);
    }
    return fluxSurfaceIntegral;
    
} 

// Evaluates the bounce average {F} of a function F = F(xi, B/Bmin) on grid point (ir,i,j). 
real_t RadialGridGenerator::CalculateBounceAverage(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType, std::function<real_t(real_t,real_t)> F){
    
    // ntheta_inter=1 signifies cylindrical geometry
    if (ntheta_interp==1){
        real_t xi0;
        if (fluxGridType == 2){
            xi0 = mg->GetXi0_f1(i,j);
        } else if (fluxGridType == 3) {
            xi0 = mg->GetXi0_f2(i,j);
        } else {
            xi0 = mg->GetXi0(i,j);
        }
        return F(xi0,1);
    } else {
    return EvaluateBounceIntegral(mg,ir,i,j,fluxGridType, F) / GetVp(mg, ir, i, j,fluxGridType);
    }
}

real_t RadialGridGenerator::EvaluateBounceIntegral(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType, std::function<real_t(real_t,real_t)> F){
    real_t 
        xi0, 
        Bmin;

    
    Bmin = this->Bmin[ir];
    xi0 = mg->GetXi0(i,j);

    if (fluxGridType == 1){
        Bmin = this->Bmin_f[ir];
    } else if (fluxGridType == 2) {
        xi0 = mg->GetXi0_f1(i,j);
    } else if (fluxGridType == 3) {
        xi0 = mg->GetXi0_f2(i,j);
    }



    std::function<real_t(real_t,real_t)> F_eff;
    
    if (GetIsTrapped(mg,ir,i,j,fluxGridType))
        F_eff = [&](real_t x, real_t  y){return ( F(x,y) + F(-x,y) )/2;};
    else 
        F_eff = F;

    real_t *B       = GetB(mg, ir,i,j,fluxGridType);
    real_t *weights = GetWeights(mg, ir,i,j,fluxGridType);
    real_t *sqrtg   = GetMetric(mg, ir,i,j,fluxGridType);

    real_t xi_particle;        
    real_t BounceIntegral = 0;
    for (len_t it = 0; it<ntheta_interp; it++) {
        xi_particle = ( (xi0 > 0) - (xi0 <0 ) ) * sqrt(1- B[it]/Bmin * (1-xi0*xi0));
        BounceIntegral += weights[it]*sqrtg[it]*F_eff(xi_particle,B[it]/Bmin);
    }        
    return BounceIntegral;
    
}

void RadialGridGenerator::EvaluateGrids(MomentumGrid **momentumGrids){
    InitializeGridQuantities(momentumGrids);

    len_t fluxGridType;
    MomentumGrid *mg;
    for (len_t ir=0; ir<nr; ir++) {
        mg = momentumGrids[ir];
        // skip to next radius if constant B, contains no trapped orbits
        // (not sure if we should be able to get here in that case)
        if (Bmin[ir]==Bmax[ir])
            continue;
        
        fluxGridType = 0;
        SetGrids(mg, ir, fluxGridType, isTrapped[ir], theta_b1[ir], theta_b2[ir], theta_bounceGrid[ir], 
        weights_bounceGrid[ir], B_bounceGrid[ir], Jacobian_bounceGrid[ir], metricSqrtG[ir], Vp[ir]);


        fluxGridType = 2;
        SetGrids(mg, ir, fluxGridType, isTrapped_f1[ir], theta_b1_f1[ir], theta_b2_f1[ir], theta_bounceGrid_f1[ir], 
        weights_bounceGrid_f1[ir], B_bounceGrid_f1[ir], Jacobian_bounceGrid_f1[ir],  metricSqrtG_f1[ir], Vp_f1[ir]);

        fluxGridType = 3;
        SetGrids(mg, ir, fluxGridType, isTrapped_f2[ir], theta_b1_f2[ir], theta_b2_f2[ir], theta_bounceGrid_f2[ir], 
        weights_bounceGrid_f2[ir], B_bounceGrid_f2[ir], Jacobian_bounceGrid_f2[ir],  metricSqrtG_f2[ir], Vp_f2[ir]);

        VpVol[ir] = 0;
        for (len_t it=0; it<ntheta_interp; it++){
            VpVol[ir] += weights[it] * Jacobian[ir][it];
        }

    }

    for (len_t ir=0; ir<nr+1; ir++) {
        
        // skip to next radius if constant B, contains no trapped orbits
        // (not sure if we should be able to get here in that case)
        if (Bmin_f[ir]==Bmax_f[ir])
            continue;
        
        fluxGridType = 1;
        SetGrids(mg, ir, fluxGridType, isTrapped_fr[ir], theta_b1_fr[ir], theta_b2_fr[ir], theta_bounceGrid_fr[ir], 
        weights_bounceGrid_fr[ir], B_bounceGrid_fr[ir], Jacobian_bounceGrid_fr[ir],  metricSqrtG_fr[ir], Vp_fr[ir]);


        VpVol_fr[ir] = 0;
        for (len_t it=0; it<ntheta_interp; it++){
            VpVol_fr[ir] += weights[it] * Jacobian_f[ir][it];
        }

    }
}



void RadialGridGenerator::SetGrids(MomentumGrid *mg, len_t ir, len_t fluxGridType, bool *isTrapped, 
    real_t *theta_b1, real_t *theta_b2, real_t **theta_bounceGrid, real_t **weights_bounceGrid, 
    real_t **B_bounceGrid, real_t **Jacobian_bounceGrid, real_t **metricSqrtG, real_t *VPrime){

    len_t np1 = mg->GetNp1();
    len_t np2 = mg->GetNp2();
    
    if (fluxGridType == 2)
        np1+=1;
    else if (fluxGridType == 3)
        np2+=1;

    real_t Bmin, Bmax;
    if (fluxGridType == 1){
        Bmin = this->Bmin_f[ir];
        Bmax = this->Bmax_f[ir];
    } else{
        Bmin = this->Bmin[ir];
        Bmax = this->Bmax[ir];
    }

    isTrapped = new bool[np1*np2];
    theta_b1  = new real_t[np1*np2];
    theta_b2  = new real_t[np1*np2];
    theta_bounceGrid    = new real_t*[np1*np2];
    weights_bounceGrid  = new real_t*[np1*np2];
    B_bounceGrid        = new real_t*[np1*np2];
    metricSqrtG = new real_t*[np1*np2];
    VPrime = new real_t[np1*np2];
    real_t xi0;
    real_t *metric_tmp;
    for (len_t i = 0; i<np1; i++){
        for (len_t j = 0; j<np2; j++){
            xi0 = mg->GetXi0(i,j);
            if (fluxGridType==2) {
                xi0 = mg->GetXi0_f1(i,j);
            } else if (fluxGridType == 3){
                xi0 = mg->GetXi0_f2(i,j);
            }

            if ( Bmax/Bmin * (1-xi0*xi0) > 1 ){
                isTrapped[j*np1+i] = true;
                SetBounceGrid(mg, ir,i,j, fluxGridType, &theta_b1[j*np1+i], &theta_b2[j*np1+i], 
                    theta_bounceGrid[j*np1+i], weights_bounceGrid[j*np1+i], B_bounceGrid[j*np1+i], 
                    Jacobian_bounceGrid[j*np1+i], metricSqrtG[j*np1+i]);
            } else {
                isTrapped[j*np1+i] = false;

                // Set metric on (passing-particle) theta grid.    
                metric_tmp = nullptr;
                mg->EvaluateMetric(i,j,fluxGridType,ntheta_interp,theta,B[ir],Bmin,metric_tmp);
                metricSqrtG[j*np1+i] = new real_t[ntheta_interp];
                for (len_t it=0; it<ntheta_interp; it++) {
                    metricSqrtG[j*np1+i][it] = Jacobian[ir][it]*metric_tmp[it];
                }

                
            }

            VPrime[j*np1+i] = EvaluateBounceIntegral(mg,ir,i,j,fluxGridType,[&](real_t,real_t){return 1;});

        }
    }

}



void RadialGridGenerator::SetBounceGrid(MomentumGrid *mg , len_t ir, len_t i, len_t j, len_t fluxGridType, real_t *theta_b1, 
                real_t *theta_b2, real_t *thetaGrid, real_t *weightsGrid, real_t *B, real_t *Jacobian, real_t *metric) {
    real_t xi0;
    if (fluxGridType == 2)
        xi0 = mg->GetXi0_f1(i,j);
    else if (fluxGridType == 3)
        xi0 = mg->GetXi0_f2(i,j);
    else 
        xi0 = mg->GetXi0(i,j);
    

    FindBouncePoints(ir,xi0,fluxGridType==1,theta_b1,theta_b2);

    gsl_spline *B_interper;
    gsl_spline *J_interper;
    if (fluxGridType==1) {
        B_interper = B_interpolator_fr[ir];
        J_interper = Jacobian_interpolator_fr[ir];
    } else {
        B_interper = B_interpolator[ir];
        J_interper = Jacobian_interpolator[ir];
    }
    thetaGrid   = new real_t[ntheta_interp];
    weightsGrid = new real_t[ntheta_interp];
    B           = new real_t[ntheta_interp];
    metric      = new real_t[ntheta_interp];

    // Linearly maps x \in [0,1] to thetaGrid \in [theta_b1, theta_b2]
    for (len_t it=0; it<ntheta_interp; it++) {
        thetaGrid[it]   = *theta_b1 + (*theta_b2-*theta_b1) * x_GL_ref[it] ;
        weightsGrid[it] = (*theta_b2 - *theta_b1) * weights_GL_ref[it];
        B[it]        = gsl_spline_eval(B_interper, thetaGrid[it], gsl_acc);
        Jacobian[it] = gsl_spline_eval(J_interper, thetaGrid[it], gsl_acc);
    }
    real_t *metric_tmp = nullptr;
    real_t Bmin;
    if (fluxGridType == 1) 
        Bmin = this->Bmin_f[ir];
    else 
        Bmin = this->Bmin[ir];

    mg->EvaluateMetric(i,j,fluxGridType,ntheta_interp,thetaGrid,B,Bmin,metric_tmp);
    metric = new real_t[ntheta_interp];
    for (len_t it=0; it<ntheta_interp; it++) {
        metric[it] = Jacobian[it]*metric_tmp[it];
    }
    

}




// Takes a theta interval theta \in [x_lower, x_upper] and iterates at most 10 times (or to a relative error of 0.001)
// and finds an estimate for the bounce point theta_bounce \in [x_lower, x_upper] (output)
void RadialGridGenerator::FindThetaBounceRoots(real_t *x_lower, real_t *x_upper, real_t *root, gsl_function gsl_func){
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    gsl_root_fsolver_set (s, &gsl_func, *x_lower, *x_upper); // finds root in [0,pi] using GSL_rootsolver_type algorithm

    int status;
    real_t rel_error = 1e-3;
    len_t max_iter = 10;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status   = gsl_root_fsolver_iterate (s);
        *root    = gsl_root_fsolver_root (s);
        *x_lower = gsl_root_fsolver_x_lower (s);
        *x_upper = gsl_root_fsolver_x_upper (s);
        status   = gsl_root_test_interval (*x_lower, *x_upper,
                                            0, rel_error);

      if (status == GSL_SUCCESS)
        gsl_root_fsolver_free(s);
        break;
    }

}

/**
 * The function is used in the evaluation of the effective passing fraction,
 * and represents x / <1-x B/Bmax>
 */
struct xiFuncParams {real_t xi0; real_t Bmin; gsl_spline *Binterper; gsl_interp_accel *gsl_acc;};
real_t RadialGridGenerator::xiParticleFunction(real_t theta, void *p){
    struct xiFuncParams *params = (struct xiFuncParams *) p;
    
    real_t xi0 = params->xi0; 
    real_t Bmin = params->Bmin;
    gsl_spline *Binterper = params->Binterper; 
    gsl_interp_accel *gsl_acc = params->gsl_acc;
    return 1 - (1-xi0*xi0) * gsl_spline_eval(Binterper, theta, gsl_acc) /Bmin ;
}

// calculates theta_bounce1 and theta_bounce2 with a root finding algorithm
void RadialGridGenerator::FindBouncePoints(len_t ir, real_t xi0, bool rFluxGrid, real_t *theta_b1, real_t *theta_b2){
    // define GSL function xi_particle as function of theta which evaluates B_interpolator(_fr) at theta.
    
    real_t Bmin;
    gsl_spline *Binterp;
    if(rFluxGrid){
        Bmin = this->Bmin_f[ir];
        Binterp = B_interpolator_fr[ir];
    } else {
        Bmin = this->Bmin[ir];
        Binterp = B_interpolator[ir];
    }
     
    xiFuncParams xi_params = {xi0,Bmin,Binterp,gsl_acc}; 
    gsl_function gsl_func;
    gsl_func.function = &(xiParticleFunction);
    gsl_func.params = &xi_params;


    // A guess can be estimated assuming that B ~ 1/[R+rcos(theta)] 
    // -- probably often a good guess. Then:
    // real_t cos_theta_b_guess = ((1-xi0*xi0)/Bmin - (1/Bmin + 1/Bmax)/2 ) * 2/( 1/Bmin - 1/Bmax );

    // Look for first root between 0 and pi
    real_t x_lower = 0.0;
    real_t x_upper = M_PI;
    real_t *root = nullptr;
    FindThetaBounceRoots(&x_lower, &x_upper, root, gsl_func);
    
    // if xi(theta = x_upper + epsilon) is real, the root 
    // corresponds to the lower bounce point theta_b1 
    if ( xiParticleFunction(x_upper,&xi_params) > 0 ){
        theta_b1 = root;
        x_lower = M_PI; 
        x_upper = 2*M_PI;
        FindThetaBounceRoots(&x_lower, &x_upper,theta_b2, gsl_func);
        //*theta_b2 = (x_lower + x_upper)/2;
    } else {
        theta_b2 = root;
        x_lower = -M_PI;
        x_upper = 0;
        FindThetaBounceRoots(&x_lower, &x_upper, theta_b1, gsl_func);
    }
}


bool RadialGridGenerator::GetIsTrapped(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType){
    len_t np1 = mg->GetNp1();

    switch (fluxGridType){
        case 0:
            return isTrapped[ir][j*np1+i];
        case 1:
            return isTrapped_fr[ir][j*np1+i];
        case 2:
            return isTrapped_f1[ir][j*(np1+1)+i];
        case 3:
            return isTrapped_f2[ir][j*np1+i];
        default:
            return NULL;
    }
}


real_t* RadialGridGenerator::GetB(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType){
    len_t np1 = mg->GetNp1();
    if (GetIsTrapped(mg,ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case 0:
                return B_bounceGrid[ir][j*np1+i];
            case 1:
                return B_bounceGrid_fr[ir][j*np1+i];
            case 2:
                return B_bounceGrid_f1[ir][j*(np1+1)+i];
            case 3:
                return B_bounceGrid_f2[ir][j*np1+i];
            default: 
                return nullptr;
        }
    } else {
        if (fluxGridType==1) {
            return B_f[ir];
        } else { 
            return B[ir];        
        }
    }
}
real_t* RadialGridGenerator::GetTheta(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType){
    len_t np1 = mg->GetNp1();
    if (GetIsTrapped(mg,ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case 0:
                return theta_bounceGrid[ir][j*np1+i];
            case 1:
                return theta_bounceGrid_fr[ir][j*np1+i];
            case 2:
                return theta_bounceGrid_f1[ir][j*(np1+1)+i];
            case 3:
                return theta_bounceGrid_f2[ir][j*np1+i];
            default: 
                return nullptr;
        }
    } else {
        return theta;
    }
}
real_t* RadialGridGenerator::GetWeights(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType){
    len_t np1 = mg->GetNp1();
    if (GetIsTrapped(mg,ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case 0:
                return weights_bounceGrid[ir][j*np1+i];
            case 1:
                return weights_bounceGrid_fr[ir][j*np1+i];
            case 2:
                return weights_bounceGrid_f1[ir][j*(np1+1)+i];
            case 3:
                return weights_bounceGrid_f2[ir][j*np1+i];
            default: 
                return nullptr;
        }
    } else {
        return weights;        
    }   
}
real_t* RadialGridGenerator::GetMetric(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType){
    len_t np1 = mg->GetNp1();
    switch (fluxGridType){
        case 0:
            return metricSqrtG[ir][j*np1+i];
        case 1:
            return metricSqrtG_fr[ir][j*np1+i];
        case 2:
            return metricSqrtG_f1[ir][j*(np1+1)+i];
        case 3:
            return metricSqrtG_f2[ir][j*np1+i];
        default: 
            return nullptr;
    }
}

real_t RadialGridGenerator::GetVp(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType){
    len_t np1 = mg->GetNp1();
    switch (fluxGridType){
        case 0:
            return Vp[ir][j*np1+i];
        case 1:
            return Vp_fr[ir][j*np1+i];
        case 2:
            return Vp_f1[ir][j*(np1+1)+i];
        case 3:
            return Vp_f2[ir][j*np1+i];
        default: 
            return -1;
    }
}

real_t *RadialGridGenerator::GetVp(len_t ir, len_t fluxGridType){
    switch (fluxGridType){
        case 0:
            return Vp[ir];
        case 1:
            return Vp_fr[ir];
        case 2:
            return Vp_f1[ir];
        case 3:
            return Vp_f2[ir];
        default: 
            return nullptr;
    }
}


real_t **RadialGridGenerator::GetVp(len_t fluxGridType){
    switch (fluxGridType){
        case 0:
            return Vp;
        case 1:
            return Vp_fr;
        case 2:
            return Vp_f1;
        case 3:
            return Vp_f2;
        default: 
            return nullptr;
    }
}





real_t RadialGridGenerator::GetVpVol(len_t ir,bool rFluxGrid){
    if (rFluxGrid)
        return VpVol_fr[ir];
    else 
        return VpVol[ir];
}

real_t *RadialGridGenerator::GetVpVol(bool rFluxGrid){
    if (rFluxGrid)
        return VpVol_fr;
    else 
        return VpVol;
}





void RadialGridGenerator::InitializeInterpolators(){
    // Below implementation uses 1D interpolation in (r,theta) data.
    // Will need to generalize to 2D interpolation somewhere.
    // will probably first implement a simple bilinear interpolation,
    // which requires B to be given in a rectilinear (r,theta) grid, i.e.
    // a tensor produced of a radial grid and a theta grid.
    // Or, just add a middle step to create B_ref and B_ref_f etc
    // in e.g. NumericalBRadialGridGenerator,
    // on which we can here do 1D interpolation? 

    for ( len_t ir = 0; nr; ir++){
        B_interpolator[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(B_interpolator[ir], theta_ref, B_ref[ir], ntheta_ref);
        Jacobian_interpolator[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(Jacobian_interpolator[ir], theta_ref, Jacobian_ref[ir], ntheta_ref);
        ROverR0_interpolator[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(ROverR0_interpolator[ir], theta_ref, ROverR0_ref[ir], ntheta_ref);
        NablaR2_interpolator[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(NablaR2_interpolator[ir], theta_ref, NablaR2_ref[ir], ntheta_ref);
    }
    for ( len_t ir = 0; nr+1; ir++){
        B_interpolator_fr[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(B_interpolator_fr[ir], theta_ref, B_ref_f[ir], ntheta_ref);
        Jacobian_interpolator_fr[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(Jacobian_interpolator_fr[ir], theta_ref, Jacobian_ref_f[ir], ntheta_ref);
        ROverR0_interpolator_fr[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(ROverR0_interpolator_fr[ir], theta_ref, ROverR0_ref_f[ir], ntheta_ref);
        NablaR2_interpolator_fr[ir] = gsl_spline_alloc(gsl_interp_steffen, ntheta_ref);
        gsl_spline_init(NablaR2_interpolator_fr[ir], theta_ref, NablaR2_ref_f[ir], ntheta_ref);
    }

    

}


void RadialGridGenerator::DeallocateInterpolators(){
    for ( len_t ir = 0; nr; ir++){
        gsl_spline_free(B_interpolator[ir]);
        gsl_spline_free(Jacobian_interpolator[ir]);
        gsl_spline_free(ROverR0_interpolator[ir]);
        gsl_spline_free(NablaR2_interpolator[ir]);
    }
    for ( len_t ir = 0; nr+1; ir++){
        gsl_spline_free(B_interpolator_fr[ir]);
        gsl_spline_free(Jacobian_interpolator_fr[ir]);
        gsl_spline_free(ROverR0_interpolator_fr[ir]);
        gsl_spline_free(NablaR2_interpolator_fr[ir]);
    }
}



void RadialGridGenerator::DeallocateMagneticFieldData(){
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
    delete [] Bmin;
    delete [] Bmin_f;
    delete [] Bmax;
    delete [] Bmax_f;
    delete [] B_ref;
    delete [] Jacobian_ref;
    delete [] ROverR0_ref;
    delete [] NablaR2_ref;
    delete [] B_ref_f;
    delete [] Jacobian_ref_f;
    delete [] ROverR0_ref_f;
    delete [] NablaR2_ref_f;
}


void RadialGridGenerator::DeallocateMagneticQuantities(){
    if (B==nullptr)
        return;

    delete [] theta;
    delete [] weights;

    for (len_t ir = 0; ir<nr; ir++){
        delete [] B[ir];
        delete [] ROverR0[ir];
        delete [] Jacobian[ir];
        delete [] NablaR2[ir];
    }
    for (len_t ir = 0; ir<nr+1; ir++){
        delete [] B_f[ir];
        delete [] ROverR0_f[ir];
        delete [] Jacobian_f[ir];
        delete [] NablaR2_f[ir];
    }

    delete [] B;
    delete [] ROverR0;
    delete [] Jacobian;
    delete [] NablaR2;
    delete [] B_f;
    delete [] ROverR0_f;
    delete [] Jacobian_f;
    delete [] NablaR2_f;

    delete [] Bmin;
    delete [] Bmax;
    delete [] Bmin_f;
    delete [] Bmax_f;
    
}
void RadialGridGenerator::InitializeMagneticQuantities(){
    DeallocateMagneticQuantities();
    
    theta   = new real_t[ntheta_interp];
    weights = new real_t[ntheta_interp];

    Bmin   = new real_t[nr];
    Bmin_f = new real_t[nr+1];
    Bmax   = new real_t[nr];
    Bmax_f = new real_t[nr+1];
    B      = new real_t*[nr];
    B_f    = new real_t*[nr+1];
    ROverR0    = new real_t*[nr];
    ROverR0_f  = new real_t*[nr+1];
    Jacobian   = new real_t*[nr];
    Jacobian_f = new real_t*[nr+1];
    NablaR2    = new real_t*[nr];
    NablaR2_f  = new real_t*[nr+1];

    
    for (len_t ir = 0; ir<nr; ir++){
        B[ir]        = new real_t[ntheta_interp];
        ROverR0[ir]  = new real_t[ntheta_interp];
        Jacobian[ir] = new real_t[ntheta_interp];
        NablaR2[ir]  = new real_t[ntheta_interp];
    }
    for (len_t ir = 0; ir<nr+1; ir++){
        B_f[ir]        = new real_t[ntheta_interp];
        ROverR0_f[ir]  = new real_t[ntheta_interp];
        Jacobian_f[ir] = new real_t[ntheta_interp];
        NablaR2_f[ir]  = new real_t[ntheta_interp];
    }
}


void RadialGridGenerator::InitializeGridQuantities(MomentumGrid **momentumGrids){
    DeallocateGridQuantities(momentumGrids);

    isTrapped    = new bool*[nr];
    isTrapped_fr = new bool*[nr+1];
    isTrapped_f1 = new bool*[nr];
    isTrapped_f2 = new bool*[nr];
    
    theta_b1    = new real_t*[nr];
    theta_b1_fr = new real_t*[nr+1];
    theta_b1_f1 = new real_t*[nr];
    theta_b1_f2 = new real_t*[nr];
    
    theta_b2    = new real_t*[nr];
    theta_b2_fr = new real_t*[nr+1];
    theta_b2_f1 = new real_t*[nr];
    theta_b2_f2 = new real_t*[nr];
    
    theta_bounceGrid    = new real_t**[nr];
    theta_bounceGrid_fr = new real_t**[nr+1];
    theta_bounceGrid_f1 = new real_t**[nr];
    theta_bounceGrid_f2 = new real_t**[nr];

    weights_bounceGrid    = new real_t**[nr];
    weights_bounceGrid_fr = new real_t**[nr+1];
    weights_bounceGrid_f1 = new real_t**[nr];
    weights_bounceGrid_f2 = new real_t**[nr];

    B_bounceGrid    = new real_t**[nr];
    B_bounceGrid_fr = new real_t**[nr+1];
    B_bounceGrid_f1 = new real_t**[nr];
    B_bounceGrid_f2 = new real_t**[nr];

    Jacobian_bounceGrid    = new real_t**[nr];
    Jacobian_bounceGrid_fr = new real_t**[nr+1];
    Jacobian_bounceGrid_f1 = new real_t**[nr];
    Jacobian_bounceGrid_f2 = new real_t**[nr];

    metricSqrtG    = new real_t**[nr];
    metricSqrtG_fr = new real_t**[nr+1];
    metricSqrtG_f1 = new real_t**[nr];
    metricSqrtG_f2 = new real_t**[nr];
    
    Vp    = new real_t*[nr];
    Vp_fr = new real_t*[nr+1];
    Vp_f1 = new real_t*[nr];
    Vp_f2 = new real_t*[nr];

    VpVol    = new real_t[nr];
    VpVol_fr = new real_t[nr+1];
}


void RadialGridGenerator::DeallocateGridQuantities(MomentumGrid **momentumGrids){
    if (isTrapped == nullptr)
        return;
    
    len_t np1, np2;
    for (len_t ir = 0; ir<nr; ir++){
        delete [] isTrapped[ir];
        delete [] isTrapped_f1[ir];
        delete [] isTrapped_f2[ir];
        /* RadialGrid owns the Vp's
        delete [] Vp[ir];
        delete [] Vp_f1[ir];
        delete [] Vp_f2[ir];
        */
        delete [] theta_b1[ir];
        delete [] theta_b1_f1[ir];
        delete [] theta_b1_f2[ir];
        delete [] theta_b2[ir];
        delete [] theta_b2_f1[ir];
        delete [] theta_b2_f2[ir];
        
        MomentumGrid *mg = momentumGrids[ir];
        np1 = mg->GetNp1();
        np2 = mg->GetNp2();
        for (len_t i = 0; i<np1; i++){
            for (len_t j = 0; j<np2; j++){
                delete [] theta_bounceGrid[ir][j*np1+i];
                delete [] weights_bounceGrid[ir][j*np1+i];
                delete [] B_bounceGrid[ir][j*np1+i];
                delete [] Jacobian_bounceGrid[ir][j*np1+i];
                delete [] metricSqrtG[ir][j*np1+i];
            }
        }
        for (len_t i = 0; i<np1+1; i++){
            for (len_t j = 0; j<np2; j++){
                delete [] theta_bounceGrid_f1[ir][j*(np1+1)+i];
                delete [] weights_bounceGrid_f1[ir][j*(np1+1)+i];
                delete [] B_bounceGrid_f1[ir][j*(np1+1)+i];
                delete [] Jacobian_bounceGrid_f1[ir][j*(np1+1)+i];
                delete [] metricSqrtG_f1[ir][j*(np1+1)+i];
            }
        }
        for (len_t i = 0; i<np1; i++){
            for (len_t j = 0; j<np2+1; j++){
                delete [] theta_bounceGrid_f2[ir][j*np1+i];
                delete [] weights_bounceGrid_f2[ir][j*np1+i];
                delete [] B_bounceGrid_f2[ir][j*np1+i];
                delete [] Jacobian_bounceGrid_f2[ir][j*np1+i];
                delete [] metricSqrtG_f2[ir][j*np1+i];
            }
        }
        delete [] theta_bounceGrid[ir];
        delete [] theta_bounceGrid_f1[ir];
        delete [] theta_bounceGrid_f2[ir];
        delete [] weights_bounceGrid[ir];
        delete [] weights_bounceGrid_f1[ir];
        delete [] weights_bounceGrid_f2[ir];
        delete [] B_bounceGrid[ir];
        delete [] B_bounceGrid_f1[ir];
        delete [] B_bounceGrid_f2[ir];
        delete [] Jacobian_bounceGrid[ir];
        delete [] Jacobian_bounceGrid_f1[ir];
        delete [] Jacobian_bounceGrid_f2[ir];
        delete [] metricSqrtG[ir];
        delete [] metricSqrtG_f1[ir];
        delete [] metricSqrtG_f2[ir];
    }

    for (len_t ir = 0; ir<nr+1; ir++){
        delete [] isTrapped_fr[ir];
        /*
        delete [] Vp_fr[ir];
        */
        delete [] theta_b1_fr[ir];
        delete [] theta_b2_fr[ir];
        
        MomentumGrid *mg = momentumGrids[ir];
        np1 = mg->GetNp1();
        np2 = mg->GetNp2();
        for (len_t i = 0; i<np1; i++){
            for (len_t j = 0; j<np2; j++){
                delete [] theta_bounceGrid_fr[ir][j*np1+i];
                delete [] weights_bounceGrid_fr[ir][j*np1+i];
                delete [] B_bounceGrid_fr[ir][j*np1+i];
                delete [] Jacobian_bounceGrid_fr[ir][j*np1+i];
                delete [] metricSqrtG_fr[ir][j*np1+i];
            }
        }
        delete [] theta_bounceGrid_fr[ir];
        delete [] weights_bounceGrid_fr[ir];
        delete [] B_bounceGrid_fr[ir];
        delete [] Jacobian_bounceGrid_fr[ir];
        delete [] metricSqrtG_fr[ir];
        
    }
    delete [] VpVol;
    delete [] VpVol_fr;

    delete [] isTrapped;
    delete [] isTrapped_fr;
    delete [] isTrapped_f1;
    delete [] isTrapped_f2;

    /* 
    delete [] Vp;
    delete [] Vp_fr;
    delete [] Vp_f1;
    delete [] Vp_f2;
    */

    delete [] theta_b1;
    delete [] theta_b1_fr;
    delete [] theta_b1_f1;
    delete [] theta_b1_f2;

    delete [] theta_b2;
    delete [] theta_b2_fr;
    delete [] theta_b2_f1;
    delete [] theta_b2_f2;
    
    delete [] theta_bounceGrid;
    delete [] theta_bounceGrid_fr;
    delete [] theta_bounceGrid_f1;
    delete [] theta_bounceGrid_f2;
    
    delete [] weights_bounceGrid;
    delete [] weights_bounceGrid_fr;
    delete [] weights_bounceGrid_f1;
    delete [] weights_bounceGrid_f2;
    
    delete [] B_bounceGrid;
    delete [] B_bounceGrid_fr;
    delete [] B_bounceGrid_f1;
    delete [] B_bounceGrid_f2;
    
    delete [] Jacobian_bounceGrid;
    delete [] Jacobian_bounceGrid_fr;
    delete [] Jacobian_bounceGrid_f1;
    delete [] Jacobian_bounceGrid_f2;
    
    delete [] metricSqrtG;
    delete [] metricSqrtG_fr;
    delete [] metricSqrtG_f1;
    delete [] metricSqrtG_f2;
    
}
