
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
    gsl_acc  = gsl_interp_accel_alloc();
}


RadialGridGenerator::~RadialGridGenerator(){
//    DeallocateGridQuantities();
    DeallocateInterpolators();
//    DeallocateMagneticFieldData();
}



void RadialGridGenerator::RebuildJacobians(RadialGrid *rGrid, MomentumGrid **momentumGrids) {
    DeallocateMagneticFieldData();
    DeallocateMagneticQuantities();
    nr = rGrid->GetNr();
    np1 = new len_t[nr+1];
    np2 = new len_t[nr+1];
    for(len_t ir=0; ir<nr; ir++){
        np1[ir] = momentumGrids[ir]->GetNp1();
        np2[ir] = momentumGrids[ir]->GetNp2();
    }
    // XXX not handling momentum grid on radial flux grid correctly
    np1[nr] = momentumGrids[0]->GetNp1(); 
    np2[nr] = momentumGrids[0]->GetNp2();

    CreateMagneticFieldData(rGrid->GetR(),rGrid->GetR_f());

    rGrid->InitializeMagneticField(
        ntheta_ref, theta_ref, B_ref, B_ref_f, Bmin, Bmin_f, Bmax, Bmax_f, Gtor, Gtor_f
    );
    
    InitializeBounceAverage(momentumGrids);

    rGrid->InitializeVprime(GetVp(FLUXGRIDTYPE_DISRIBUTION),GetVp(FLUXGRIDTYPE_RADIAL),
                            GetVp(FLUXGRIDTYPE_P1),GetVp(FLUXGRIDTYPE_P2),
                            GetVpVol(false), GetVpVol(true));


    
}

void RadialGridGenerator::InitializeBounceAverage(MomentumGrid **momentumGrids){
    if(ntheta_ref==1){
        ntheta_interp = 1;
        x_GL_ref = new real_t[1];
        weights_GL_ref = new real_t[1];

        x_GL_ref[0] = 0;
        weights_GL_ref[0] = 2*M_PI;
    } else {
        // Create reference Gauss-Legendre quadrature on x\in[0,1]
        gsl_integration_fixed_workspace *gsl_GL = gsl_integration_fixed_alloc(thetaGridType,ntheta_interp,0,1,0,0);
        this->x_GL_ref = gsl_GL->x;
        this->weights_GL_ref = gsl_GL->weights;

        InitializeInterpolators();

    } 
    InitializeMagneticQuantities();


    
    // if ntheta_ref = 1, we assume cylindrical geometry and will set
    // bounce and flux surface averaging to be identity functions
    // and avoid defining a bunch of stuff.
    if(ntheta_ref==1){
        theta[0]   = 0;
        weights[0] = 1;
        for (len_t ir=0;ir<nr;ir++){
            B[ir][0] = B_ref[ir][0];
            ROverR0[ir][0]  = ROverR0_ref[ir][0];
            Jacobian[ir][0] = Jacobian_ref[ir][0];
            NablaR2[ir][0]  = NablaR2_ref[ir][0];
            
        }
        for (len_t ir=0;ir<nr+1;ir++){
            B_f[ir][0] = B_ref_f[ir][0];
            ROverR0_f[ir][0]  = ROverR0_ref_f[ir][0];
            Jacobian_f[ir][0] = Jacobian_ref_f[ir][0];
            NablaR2_f[ir][0]  = NablaR2_ref_f[ir][0];
        }

    } else {
        real_t theta_max;
        if(isUpDownSymmetric)
            theta_max = M_PI;
        else 
            theta_max = 2*M_PI;

        // Weights are still 2pi since we multiply integral by 2 in case of updownsymmetric
        for (len_t it=0; it<ntheta_interp; it++) {
            theta[it]   = theta_max * x_GL_ref[it];
            weights[it] = 2*M_PI * weights_GL_ref[it];
        }
        
        for (len_t ir=0; ir<nr;ir++){
            for (len_t it=0; it<ntheta_interp; it++){
                B[ir][it]        = gsl_spline_eval(B_interpolator[ir], theta[it], gsl_acc);
                ROverR0[ir][it]  = gsl_spline_eval(ROverR0_interpolator[ir], theta[it], gsl_acc);
                Jacobian[ir][it] = gsl_spline_eval(Jacobian_interpolator[ir], theta[it], gsl_acc);
                NablaR2[ir][it]  = gsl_spline_eval(NablaR2_interpolator[ir], theta[it], gsl_acc);
            }
        }

        for (len_t ir=0; ir<nr+1;ir++){
            for (len_t it=0; it<ntheta_interp; it++){
                B_f[ir][it]        = gsl_spline_eval(B_interpolator_fr[ir], theta[it], gsl_acc);
                ROverR0_f[ir][it]  = gsl_spline_eval(ROverR0_interpolator_fr[ir], theta[it], gsl_acc);
                NablaR2_f[ir][it]  = gsl_spline_eval(NablaR2_interpolator_fr[ir], theta[it], gsl_acc);
                Jacobian_f[ir][it] = gsl_spline_eval(Jacobian_interpolator_fr[ir], theta[it], gsl_acc);
            }
        }
    }
    CalculateQuantities(momentumGrids);

//    DeallocateInterpolators();
}





// Evaluates the flux surface average <F> of a function F = F(B/Bmin, R/R0, |nabla r|^2) on radial grid point ir. 
real_t RadialGridGenerator::CalculateFluxSurfaceAverage(len_t ir, bool rFluxGrid, std::function<real_t(real_t,real_t,real_t)> F){
    if (ntheta_interp == 1){
        return F(1.0,1.0,1.0);
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
        fluxSurfaceIntegral += 2*M_PI*weights[it] * Jacobian[it] * F(B[it]/Bmin, ROverR0[it], NablaR2[it]);
    }
    return fluxSurfaceIntegral;
    
} 

// Evaluates the bounce average {F} of a function F = F(xi/xi0, B/Bmin) on grid point (ir,i,j). 
real_t RadialGridGenerator::CalculateBounceAverage(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t)> F){
    // ntheta_inter=1 signifies cylindrical geometry
    if (ntheta_interp==1){
        return F(1,1);
    } 

    return EvaluateBounceIntegral(mg,ir,i,j,fluxGridType, F) / GetVp(mg, ir, i, j,fluxGridType);
}



real_t RadialGridGenerator::EvaluateBounceIntegral(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t)> F){
    real_t 
        xi0, 
        Bmin;

    if (fluxGridType == FLUXGRIDTYPE_RADIAL)
        Bmin = this->Bmin_f[ir];
    else
        Bmin = this->Bmin[ir];
    

    if (fluxGridType == FLUXGRIDTYPE_P1) {
        xi0 = mg->GetXi0_f1(i,j);
    } else if (fluxGridType == FLUXGRIDTYPE_P2) {
        xi0 = mg->GetXi0_f2(i,j);
    } else {
        xi0 = mg->GetXi0(i,j);
    }
    
    std::function<real_t(real_t,real_t)> F_eff;
    // If trapped, adds contribution from -xi0, since negative xi0 are presumably not kept on the grid.
    if (GetIsTrapped(mg,ir,i,j,fluxGridType))
        F_eff = [&](real_t x, real_t  y){return  F(x,y) + F(-x,y) ;};
    else 
        F_eff = F;

    real_t *B       = GetB(mg, ir,i,j,fluxGridType);
    real_t *weights = GetWeights(mg, ir,i,j,fluxGridType);
    real_t *sqrtg   = GetMetric(mg, ir,i,j,fluxGridType);

    if(ntheta_interp==1)
        return 2*M_PI * sqrtg[0] * F_eff(1,1);

    real_t xiOverXi0,w,BOverBmin;        
    real_t BounceIntegral = 0;
    for (len_t it = 0; it<ntheta_interp; it++) {
        real_t xi2 = 1- B[it]/Bmin * (1-xi0*xi0);
        if(xi2==0){
            // TODO: this strictly only works for up-down symmetric magnetic fields, but I anticipate
            // that this will never cause unphysical behavior at xi0=0 even for asymmetric fields
            real_t thetaOverThetaBounce = it/(ntheta_interp-1);   
            xiOverXi0 = sqrt(1-thetaOverThetaBounce*thetaOverThetaBounce);
            BOverBmin = 1;
            w = 1.0/ntheta_interp;
        } else {
            xiOverXi0 = sqrt(xi2/(xi0*xi0));
            w = weights[it];
            BOverBmin = B[it]/Bmin;
        }

        BounceIntegral += 2*M_PI*w*sqrtg[it]*F_eff(xiOverXi0,BOverBmin);
    }        
    return BounceIntegral;
    
}



real_t RadialGridGenerator::evaluateXiAtTheta(len_t ir, real_t xi0, real_t theta, bool rFluxGrid){
    real_t sgnXi = (xi0>0) - (xi0<0);
    real_t Bmin;
    if(rFluxGrid){
        Bmin = this->Bmin_f[ir];
    } else{ 
        Bmin = this->Bmin[ir];
    }
    real_t BOverBmin = evaluateBAtTheta(ir,theta,rFluxGrid)/Bmin;
    return sgnXi*sqrt(1-BOverBmin*(1-xi0*xi0));
}

struct generalBounceIntegralParams {len_t ir; real_t xi0; real_t p; bool rFluxGrid; real_t Bmin; std::function<real_t(real_t,real_t)> F_eff; RadialGridGenerator *rgg; /*MomentumGrid *mg;*/};
real_t generalBounceIntegralFunc(real_t theta, void *par){
    struct generalBounceIntegralParams *params = (struct generalBounceIntegralParams *) par;
    
    len_t ir = params->ir;
    real_t xi0 = params->xi0;
    real_t p = params->p;
    bool rFluxGrid = params->rFluxGrid;
    real_t Bmin = params->Bmin;
    RadialGridGenerator *rgg = params->rgg; 
//    MomentumGrid *mg = params->mg;
    std::function<real_t(real_t,real_t)> F_eff = params->F_eff;
    real_t B = rgg->evaluateBAtTheta(ir,theta,rFluxGrid);
    real_t Jacobian = rgg->evaluateJacobianAtTheta(ir,theta,rFluxGrid);
    real_t sqrtG = MomentumGrid::evaluatePXiMetricOverP2(p,xi0,B,Bmin);
    
    real_t F =  F_eff(rgg->evaluateXiAtTheta(ir,xi0,theta,rFluxGrid)/xi0,B/Bmin);
    return 2*M_PI*Jacobian*sqrtG*F;
}


/**
 * Returns lim_{p\to 0} Vp(r,p,xi0). rFluxGrid specifies
 */
real_t RadialGridGenerator::evaluatePXiVpOverP2(len_t ir, real_t xi0, bool rFluxGrid,gsl_integration_workspace *gsl_ad_w){
    std::function<real_t(real_t,real_t)> FUnity = [](real_t,real_t){return 1;};
    real_t p = 0;
    return evaluatePXiBounceIntegralAtP(ir,  p,  xi0,  rFluxGrid, FUnity, gsl_ad_w);
}

// Evaluates the bounce average of the function F = F(xi/xi0, B/Bmin) at  
// radial grid point ir, momentum p and pitch xi0, using an adaptive quadrature.
real_t RadialGridGenerator::evaluatePXiBounceIntegralAtP(len_t ir, real_t p, real_t xi0, bool rFluxGrid, std::function<real_t(real_t,real_t)> F,gsl_integration_workspace *gsl_ad_w){
    real_t Bmin,Bmax;
    if(rFluxGrid){
        Bmin = this->Bmin_f[ir];
        Bmax = this->Bmax_f[ir];
    } else{ 
        Bmin = this->Bmin[ir];
        Bmax = this->Bmax[ir];
    }

    std::function<real_t(real_t,real_t)> F_eff;
    bool isTrapped = (Bmax/Bmin * (1-xi0*xi0) > 1);
    // If trapped, adds contribution from -xi0, since negative xi0 are presumably not kept on the grid.
    real_t theta_b1, theta_b2;
    if (isTrapped){
        F_eff = [&](real_t x, real_t  y){return  F(x,y) + F(-x,y) ;};
        FindBouncePoints(ir, xi0, rFluxGrid, &theta_b1, &theta_b2);
    } else { 
        F_eff = F;
        theta_b1 = 0;
        theta_b2 = 2*M_PI;
    }

    gsl_function GSL_func;
    GSL_func.function = &(generalBounceIntegralFunc);
    generalBounceIntegralParams params = {ir,xi0,p,rFluxGrid,Bmin,F_eff,this};
    GSL_func.params = &params;
    real_t bounceIntegral, error; 

    real_t epsabs = 0, epsrel = 1e-4, lim = gsl_ad_w->limit; 
    gsl_integration_qags(&GSL_func,theta_b1,theta_b2,epsabs,epsrel,lim,gsl_ad_w,&bounceIntegral,&error);
    return bounceIntegral;
}
// Evaluates the bounce average {F} of a function F = F(xi/xi0, B/Bmin) on grid point (ir,i,j). 
real_t RadialGridGenerator::evaluatePXiBounceAverageAtP(len_t ir, real_t p, real_t xi0, bool rFluxGrid, std::function<real_t(real_t,real_t)> F,gsl_integration_workspace *gsl_ad_w){
    
    // ntheta_interp=1 signifies cylindrical geometry
    if (ntheta_interp==1){
        return F(1,1);
    } else {
        std::function<real_t(real_t,real_t)> FUnity = [](real_t, real_t){
        return 1;};
        return evaluatePXiBounceIntegralAtP(ir, p, xi0, rFluxGrid, F,gsl_ad_w) / evaluatePXiBounceIntegralAtP(ir, p, xi0, rFluxGrid, FUnity,gsl_ad_w);
    }
}

void RadialGridGenerator::CalculateQuantities(MomentumGrid **momentumGrids){
    InitializeGridQuantities();

    fluxGridType fluxGridType;
    MomentumGrid *mg;
    for (len_t ir=0; ir<nr; ir++) {
        mg = momentumGrids[ir];
        
        fluxGridType = FLUXGRIDTYPE_DISRIBUTION;
        SetQuantities(mg, ir, fluxGridType, isTrapped, theta_b1, theta_b2, theta_bounceGrid, 
        weights_bounceGrid, B_bounceGrid, B, Jacobian, Jacobian_bounceGrid, metricSqrtG, Vp);

        fluxGridType = FLUXGRIDTYPE_P1;
        SetQuantities(mg, ir, fluxGridType, isTrapped_f1, theta_b1_f1, theta_b2_f1, theta_bounceGrid_f1, 
        weights_bounceGrid_f1, B_bounceGrid_f1, B, Jacobian, Jacobian_bounceGrid_f1,  metricSqrtG_f1, Vp_f1);

        fluxGridType = FLUXGRIDTYPE_P2;
        SetQuantities(mg, ir, fluxGridType, isTrapped_f2, theta_b1_f2, theta_b2_f2, theta_bounceGrid_f2, 
        weights_bounceGrid_f2, B_bounceGrid_f2, B, Jacobian, Jacobian_bounceGrid_f2,  metricSqrtG_f2, Vp_f2);

        VpVol[ir] = EvaluateFluxSurfaceIntegral(ir,false,[](real_t,real_t,real_t ){return 1;});
    }

    /*
    XXX I explicitly assume that momentumGrids are the same on all radii in the below block
    */ 
    for (len_t ir=0; ir<nr+1; ir++) {
        mg = momentumGrids[0];

        fluxGridType = FLUXGRIDTYPE_RADIAL;
        SetQuantities(mg, ir, fluxGridType, isTrapped_fr, theta_b1_fr, theta_b2_fr, theta_bounceGrid_fr, 
        weights_bounceGrid_fr, B_bounceGrid_fr, B_f, Jacobian_f, Jacobian_bounceGrid_fr,  metricSqrtG_fr, Vp_fr);

        VpVol_fr[ir] = EvaluateFluxSurfaceIntegral(ir,true,[](real_t,real_t,real_t ){return 1;});
    }
}



void RadialGridGenerator::SetQuantities(MomentumGrid *mg, len_t ir, fluxGridType fluxGridType, bool **&isTrapped, 
    real_t **&theta_b1, real_t **&theta_b2, real_t ***&theta_bounceGrid, real_t ***&weights_bounceGrid, 
    real_t ***&B_bounceGrid, real_t **&B, real_t **&Jacobian, real_t ***&Jacobian_bounceGrid, real_t ***&metricSqrtG, real_t **&VPrime){

    len_t np1 = mg->GetNp1();
    len_t np2 = mg->GetNp2();
    
    if (fluxGridType == FLUXGRIDTYPE_P1)
        np1+=1;
    else if (fluxGridType == FLUXGRIDTYPE_P2)
        np2+=1;

    real_t Bmin, Bmax;
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        Bmin = this->Bmin_f[ir];
        Bmax = this->Bmax_f[ir];
    } else{
        Bmin = this->Bmin[ir];
        Bmax = this->Bmax[ir];
    }

    isTrapped[ir] = new bool[np1*np2];
    theta_b1[ir]  = new real_t[np1*np2];
    theta_b2[ir]  = new real_t[np1*np2];
    theta_bounceGrid[ir]    = new real_t*[np1*np2];
    weights_bounceGrid[ir]  = new real_t*[np1*np2];
    B_bounceGrid[ir]        = new real_t*[np1*np2];
    Jacobian_bounceGrid[ir] = new real_t*[np1*np2];
    metricSqrtG[ir] = new real_t*[np1*np2];
    VPrime[ir] = new real_t[np1*np2];
    std::function<real_t(real_t,real_t)> FUnity = [](real_t,real_t){return 1;};
    real_t xi0;
    len_t ind;
    for (len_t i = 0; i<np1; i++){
        for (len_t j = 0; j<np2; j++){
            
            if (fluxGridType==FLUXGRIDTYPE_P1) {
                xi0 = mg->GetXi0_f1(i,j);
            } else if (fluxGridType == FLUXGRIDTYPE_P2){
                xi0 = mg->GetXi0_f2(i,j);
            } else {
                xi0 = mg->GetXi0(i,j);
            }
            ind = j*np1+i;
            metricSqrtG[ir][ind] = new real_t[ntheta_interp];
            if ( Bmax/Bmin * (1-xi0*xi0) > 1 ){
                isTrapped[ir][ind] = true;
                SetBounceGrid(mg, ir,i,j, fluxGridType, theta_b1, theta_b2, 
                    theta_bounceGrid, weights_bounceGrid, B_bounceGrid, 
                    Jacobian_bounceGrid, metricSqrtG);
            } else {
                isTrapped[ir][ind] = false;
                // Set metric on (passing-particle) theta grid.    
                
                mg->EvaluateMetric(i,j,fluxGridType,ntheta_interp,theta,B[ir],Bmin,metricSqrtG[ir][ind]);
                for (len_t it=0; it<ntheta_interp; it++) {
                    metricSqrtG[ir][ind][it] *= Jacobian[ir][it];
                }
            }

            VPrime[ir][ind] = EvaluateBounceIntegral(mg,ir,i,j,fluxGridType,FUnity);
        }
    }

}


void RadialGridGenerator::SetBounceGrid(MomentumGrid *mg , len_t ir, len_t i, len_t j, fluxGridType fluxGridType, real_t **&theta_b1, 
                real_t **&theta_b2, real_t ***&thetaGrid, real_t ***&weightsGrid, real_t ***&B, real_t ***&Jacobian, real_t ***&metric) {
    real_t xi0;
    if (fluxGridType == FLUXGRIDTYPE_P1)
        xi0 = mg->GetXi0_f1(i,j);
    else if (fluxGridType == FLUXGRIDTYPE_P2)
        xi0 = mg->GetXi0_f2(i,j);
    else 
        xi0 = mg->GetXi0(i,j);

    len_t ind = j*(mg->GetNp1()+(fluxGridType==FLUXGRIDTYPE_P1))+i; 
    FindBouncePoints(ir,xi0,fluxGridType==FLUXGRIDTYPE_RADIAL,&theta_b1[ir][ind],&theta_b2[ir][ind]);

    gsl_spline *B_interper;
    gsl_spline *J_interper;
    if (fluxGridType==FLUXGRIDTYPE_RADIAL) {
        B_interper = B_interpolator_fr[ir];
        J_interper = Jacobian_interpolator_fr[ir];
    } else {
        B_interper = B_interpolator[ir];
        J_interper = Jacobian_interpolator[ir];
    }
    thetaGrid[ir][ind]   = new real_t[ntheta_interp];
    weightsGrid[ir][ind] = new real_t[ntheta_interp];
    B[ir][ind]           = new real_t[ntheta_interp];
    Jacobian[ir][ind]    = new real_t[ntheta_interp];
    real_t *metric_tmp   = new real_t[ntheta_interp];

    // if symmetric flux surface, take grid from 0 to upper bounce point theta_b2, and 
    // multiply quadrature weights by 2 
    if (isUpDownSymmetric){
        for (len_t it=0; it<ntheta_interp; it++) {
            thetaGrid[ir][ind][it]   = theta_b2[ir][ind] * x_GL_ref[it] ;
            weightsGrid[ir][ind][it] = 2 * theta_b2[ir][ind] * weights_GL_ref[it];
            B[ir][ind][it]        = gsl_spline_eval(B_interper, thetaGrid[ir][ind][it], gsl_acc);
            Jacobian[ir][ind][it] = gsl_spline_eval(J_interper, thetaGrid[ir][ind][it], gsl_acc);
        }
    } else{
    // Linearly maps x \in [0,1] to thetaGrid \in [theta_b1, theta_b2]
        for (len_t it=0; it<ntheta_interp; it++) {
            thetaGrid[ir][ind][it]   = theta_b1[ir][ind] + (theta_b2[ir][ind]-theta_b1[ir][ind]) * x_GL_ref[it] ;
            weightsGrid[ir][ind][it] = (theta_b2[ir][ind] - theta_b1[ir][ind]) * weights_GL_ref[it];
            B[ir][ind][it]        = gsl_spline_eval(B_interper, thetaGrid[ir][ind][it], gsl_acc);
            Jacobian[ir][ind][it] = gsl_spline_eval(J_interper, thetaGrid[ir][ind][it], gsl_acc);
        }
    }
    real_t Bmin;
    if (fluxGridType == FLUXGRIDTYPE_RADIAL) 
        Bmin = this->Bmin_f[ir];
    else 
        Bmin = this->Bmin[ir];
    mg->EvaluateMetric(i,j,fluxGridType,ntheta_interp,thetaGrid[ir][ind],B[ir][ind],Bmin,metric_tmp);
    metric[ir][ind] = metric_tmp;
    for (len_t it=0; it<ntheta_interp; it++) {
        metric[ir][ind][it] *= Jacobian[ir][ind][it];
    }

}




// Takes a theta interval theta \in [x_lower, x_upper] and iterates at most max_iter (=10) times (or to a relative error of 0.001)
// to find an estimate for the bounce point theta_bounce \in [x_lower, x_upper]
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

        if (status == GSL_SUCCESS){
            gsl_root_fsolver_free(s);
            break;
        }
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
    real_t root = 0;
    FindThetaBounceRoots(&x_lower, &x_upper, &root, gsl_func);
    
    // In symmetric field, this root corresponds to theta_b2
    if (isUpDownSymmetric){
        *theta_b2 = root;
        *theta_b1 = -root;
    } else {

        // if xi(theta = x_upper + epsilon) is real, the root 
        // corresponds to the lower bounce point theta_b1 
        if ( xiParticleFunction(x_upper,&xi_params) > 0 ){
            theta_b1 = &root;
            x_lower = M_PI; 
            x_upper = 2*M_PI;
            FindThetaBounceRoots(&x_lower, &x_upper,theta_b2, gsl_func);
            //*theta_b2 = (x_lower + x_upper)/2;
        } else {
            theta_b2 = &root;
            x_lower = -M_PI;
            x_upper = 0;
            FindThetaBounceRoots(&x_lower, &x_upper, theta_b1, gsl_func);
        }
    }
}


bool RadialGridGenerator::GetIsTrapped(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType){
    len_t np1 = mg->GetNp1();
    len_t ind = j*np1+i;
    switch (fluxGridType){
        case FLUXGRIDTYPE_DISRIBUTION:
            return isTrapped[ir][ind];
        case FLUXGRIDTYPE_RADIAL:
            return isTrapped_fr[ir][ind];
        case FLUXGRIDTYPE_P1:
            return isTrapped_f1[ir][ind+j];
        case FLUXGRIDTYPE_P2:
            return isTrapped_f2[ir][ind];
        default:
            return NULL;
    }
}


real_t* RadialGridGenerator::GetB(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType){
    len_t np1 = mg->GetNp1();
    len_t ind = j*np1+i;
    
    if (GetIsTrapped(mg,ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISRIBUTION:
                return B_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return B_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return B_bounceGrid_f1[ir][ind+j];
            case FLUXGRIDTYPE_P2:
                return B_bounceGrid_f2[ir][ind];
            default: 
                return nullptr;
        }
    } else {
        if (fluxGridType==FLUXGRIDTYPE_RADIAL) {
            return B_f[ir];
        } else { 
            return B[ir];        
        }
    }
}
real_t* RadialGridGenerator::GetTheta(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType){
    len_t np1 = mg->GetNp1();
    len_t ind = j*np1+i;
    
    if (GetIsTrapped(mg,ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISRIBUTION:
                return theta_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return theta_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return theta_bounceGrid_f1[ir][ind+j];
            case FLUXGRIDTYPE_P2:
                return theta_bounceGrid_f2[ir][ind];
            default: 
                return nullptr;
        }
    } else {
        return theta;
    }
}
real_t* RadialGridGenerator::GetWeights(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType){
    len_t np1 = mg->GetNp1();
    len_t ind = j*np1+i;
    if (GetIsTrapped(mg,ir,i,j,fluxGridType)) {
        switch (fluxGridType){
            case FLUXGRIDTYPE_DISRIBUTION:
                return weights_bounceGrid[ir][ind];
            case FLUXGRIDTYPE_RADIAL:
                return weights_bounceGrid_fr[ir][ind];
            case FLUXGRIDTYPE_P1:
                return weights_bounceGrid_f1[ir][ind+j];
            case FLUXGRIDTYPE_P2:
                return weights_bounceGrid_f2[ir][ind];
            default: 
                return nullptr;
        }
    } else {
        return weights;        
    }   
}
real_t* RadialGridGenerator::GetMetric(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType){
    len_t np1 = mg->GetNp1();
    len_t ind = j*np1+i;
    switch (fluxGridType){
        case FLUXGRIDTYPE_DISRIBUTION:
            return metricSqrtG[ir][ind];
        case FLUXGRIDTYPE_RADIAL:
            return metricSqrtG_fr[ir][ind];
        case FLUXGRIDTYPE_P1:
            return metricSqrtG_f1[ir][ind+j];
        case FLUXGRIDTYPE_P2:
            return metricSqrtG_f2[ir][ind];
        default: 
            return nullptr;
    }
}

real_t RadialGridGenerator::GetVp(MomentumGrid *mg, len_t ir, len_t i, len_t j, fluxGridType fluxGridType){
    len_t np1 = mg->GetNp1();
    len_t ind = j*np1+i;
    switch (fluxGridType){
        case FLUXGRIDTYPE_DISRIBUTION:
            return Vp[ir][ind];
        case FLUXGRIDTYPE_RADIAL:
            return Vp_fr[ir][ind];
        case FLUXGRIDTYPE_P1:
            return Vp_f1[ir][ind+j];
        case FLUXGRIDTYPE_P2:
            return Vp_f2[ir][ind];
        default: 
            return -1;
    }
}

real_t *RadialGridGenerator::GetVp(len_t ir, fluxGridType fluxGridType){
    switch (fluxGridType){
        case FLUXGRIDTYPE_DISRIBUTION:
            return Vp[ir];
        case FLUXGRIDTYPE_RADIAL:
            return Vp_fr[ir];
        case FLUXGRIDTYPE_P1:
            return Vp_f1[ir];
        case FLUXGRIDTYPE_P2:
            return Vp_f2[ir];
        default: 
            return nullptr;
    }
}


real_t **RadialGridGenerator::GetVp(fluxGridType fluxGridType){
    switch (fluxGridType){
        case FLUXGRIDTYPE_DISRIBUTION:
            return Vp;
        case FLUXGRIDTYPE_RADIAL:
            return Vp_fr;
        case FLUXGRIDTYPE_P1:
            return Vp_f1;
        case FLUXGRIDTYPE_P2:
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

real_t RadialGridGenerator::evaluateBAtTheta(len_t ir, real_t theta, bool rFluxGrid){
    if (theta < 0)
        theta += 2*M_PI;
    if(rFluxGrid){
        return gsl_spline_eval(B_interpolator_fr[ir], theta, gsl_acc);
    } else {
        return gsl_spline_eval(B_interpolator[ir], theta, gsl_acc);
    }
}

real_t RadialGridGenerator::evaluateJacobianAtTheta(len_t ir, real_t theta, bool rFluxGrid){
    if (theta < 0)
        theta += 2*M_PI;
    if(rFluxGrid){
        return gsl_spline_eval(Jacobian_interpolator_fr[ir], theta, gsl_acc);
    } else {
        return gsl_spline_eval(Jacobian_interpolator[ir], theta, gsl_acc);
    }
}

real_t RadialGridGenerator::evaluateROverR0AtTheta(len_t ir, real_t theta, bool rFluxGrid){
    if (theta < 0)
        theta += 2*M_PI;
    if(rFluxGrid){
        return gsl_spline_eval(ROverR0_interpolator_fr[ir], theta, gsl_acc);
    } else {
        return gsl_spline_eval(ROverR0_interpolator[ir], theta, gsl_acc);
    }
}

real_t RadialGridGenerator::evaluateNablaR2AtTheta(len_t ir, real_t theta, bool rFluxGrid){
    if (theta < 0)
        theta += 2*M_PI;
    if(rFluxGrid){
        return gsl_spline_eval(NablaR2_interpolator_fr[ir], theta, gsl_acc);
    } else {
        return gsl_spline_eval(NablaR2_interpolator[ir], theta, gsl_acc);
    }
}





void RadialGridGenerator::InitializeInterpolators(){
    DeallocateInterpolators();
    // Below implementation uses 1D interpolation in (r,theta) data.
    // Will need to generalize to 2D interpolation somewhere.
    // will probably first implement a simple bilinear interpolation,
    // which requires B to be given in a rectilinear (r,theta) grid, i.e.
    // a tensor produced of a radial grid and a theta grid.
    // Or, just add a middle step to create B_ref and B_ref_f etc
    // in e.g. NumericalBRadialGridGenerator,
    // on which we can here do 1D interpolation? 

    B_interpolator           = new gsl_spline*[nr];
    B_interpolator_fr        = new gsl_spline*[nr+1];
    Jacobian_interpolator    = new gsl_spline*[nr];
    Jacobian_interpolator_fr = new gsl_spline*[nr+1];
    ROverR0_interpolator     = new gsl_spline*[nr];
    ROverR0_interpolator_fr  = new gsl_spline*[nr+1];
    NablaR2_interpolator     = new gsl_spline*[nr];
    NablaR2_interpolator_fr  = new gsl_spline*[nr+1];
    const gsl_interp_type *gsl_type = gsl_interp_linear;
//    const gsl_interp_type *gsl_type = gsl_interp_steffen;
    for ( len_t ir = 0; ir<nr; ir++){
        B_interpolator[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(B_interpolator[ir], theta_ref, B_ref[ir], ntheta_ref);
        Jacobian_interpolator[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(Jacobian_interpolator[ir], theta_ref, Jacobian_ref[ir], ntheta_ref);
        ROverR0_interpolator[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(ROverR0_interpolator[ir], theta_ref, ROverR0_ref[ir], ntheta_ref);
        NablaR2_interpolator[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(NablaR2_interpolator[ir], theta_ref, NablaR2_ref[ir], ntheta_ref);
    }
    for ( len_t ir = 0; ir<nr+1; ir++){
        B_interpolator_fr[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(B_interpolator_fr[ir], theta_ref, B_ref_f[ir], ntheta_ref);
        Jacobian_interpolator_fr[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(Jacobian_interpolator_fr[ir], theta_ref, Jacobian_ref_f[ir], ntheta_ref);
        ROverR0_interpolator_fr[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(ROverR0_interpolator_fr[ir], theta_ref, ROverR0_ref_f[ir], ntheta_ref);
        NablaR2_interpolator_fr[ir] = gsl_spline_alloc(gsl_type, ntheta_ref);
        gsl_spline_init(NablaR2_interpolator_fr[ir], theta_ref, NablaR2_ref_f[ir], ntheta_ref);
    }

    

}


void RadialGridGenerator::DeallocateInterpolators(){
    if (B_interpolator == nullptr)
        return;
    for ( len_t ir = 0; ir<nr; ir++){
        gsl_spline_free(B_interpolator[ir]);
        gsl_spline_free(Jacobian_interpolator[ir]);
        gsl_spline_free(ROverR0_interpolator[ir]);
        gsl_spline_free(NablaR2_interpolator[ir]);
    }
    for ( len_t ir = 0; ir<nr+1; ir++){
        gsl_spline_free(B_interpolator_fr[ir]);
        gsl_spline_free(Jacobian_interpolator_fr[ir]);
        gsl_spline_free(ROverR0_interpolator_fr[ir]);
        gsl_spline_free(NablaR2_interpolator_fr[ir]);
    }
    delete [] B_interpolator;
    delete [] B_interpolator_fr;
    delete [] Jacobian_interpolator;
    delete [] Jacobian_interpolator_fr;
    delete [] ROverR0_interpolator;
    delete [] ROverR0_interpolator_fr;
    delete [] NablaR2_interpolator;
    delete [] NablaR2_interpolator_fr;

}



void RadialGridGenerator::DeallocateMagneticFieldData(){
    if(np1 != nullptr){
        delete [] np1;
        delete [] np2;
    }
    if (B_ref==nullptr)
        return;

    for(len_t ir = 0; ir<nr; ir++){
        delete [] B_ref[ir];
        delete [] Jacobian_ref[ir];
        delete [] ROverR0_ref[ir];
        delete [] NablaR2_ref[ir];
    }
    for(len_t ir = 0; ir<nr; ir++){
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
//    delete [] Gtor;
//    delete [] Gtor_f;
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
/*
    delete [] Bmin;
    delete [] Bmax;
    delete [] Bmin_f;
    delete [] Bmax_f;
*/  
}
void RadialGridGenerator::InitializeMagneticQuantities(){
    DeallocateMagneticQuantities();
    
    theta   = new real_t[ntheta_interp];
    weights = new real_t[ntheta_interp];

/*
    Bmin   = new real_t[nr];
    Bmin_f = new real_t[nr+1];
    Bmax   = new real_t[nr];
    Bmax_f = new real_t[nr+1];
*/
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


void RadialGridGenerator::InitializeGridQuantities(){
    DeallocateGridQuantities();

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


void RadialGridGenerator::DeallocateGridQuantities(){
    if (isTrapped == nullptr)
        return;
    
    len_t ind;
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
        
        for (len_t i = 0; i<np1[ir]; i++){
            for (len_t j = 0; j<np2[ir]; j++){
                ind = j*np1[ir]+i;
                delete [] theta_bounceGrid[ir][ind];
                delete [] weights_bounceGrid[ir][ind];
                delete [] B_bounceGrid[ir][ind];
                delete [] Jacobian_bounceGrid[ir][ind];
                delete [] metricSqrtG[ir][ind];
            }
        }
        for (len_t i = 0; i<np1[ir]+1; i++){
            for (len_t j = 0; j<np2[ir]; j++){
                ind = j*(np1[ir]+1)+i;
                delete [] theta_bounceGrid_f1[ir][ind];
                delete [] weights_bounceGrid_f1[ir][ind];
                delete [] B_bounceGrid_f1[ir][ind];
                delete [] Jacobian_bounceGrid_f1[ir][ind];
                delete [] metricSqrtG_f1[ir][ind];
            }
        }
        for (len_t i = 0; i<np1[ir]; i++){
            for (len_t j = 0; j<np2[ir]+1; j++){
                ind = j*np1[ir]+i;
                delete [] theta_bounceGrid_f2[ir][ind];
                delete [] weights_bounceGrid_f2[ir][ind];
                delete [] B_bounceGrid_f2[ir][ind];
                delete [] Jacobian_bounceGrid_f2[ir][ind];
                delete [] metricSqrtG_f2[ir][ind];
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
        
        for (len_t i = 0; i<np1[ir]; i++){
            for (len_t j = 0; j<np2[ir]; j++){
                ind = j*np1[ir]+i;
                delete [] theta_bounceGrid_fr[ir][ind];
                delete [] weights_bounceGrid_fr[ir][ind];
                delete [] B_bounceGrid_fr[ir][ind];
                delete [] Jacobian_bounceGrid_fr[ir][ind];
                delete [] metricSqrtG_fr[ir][ind];
            }
        }
        delete [] theta_bounceGrid_fr[ir];
        delete [] weights_bounceGrid_fr[ir];
        delete [] B_bounceGrid_fr[ir];
        delete [] Jacobian_bounceGrid_fr[ir];
        delete [] metricSqrtG_fr[ir];
        
    }
    /*
    delete [] VpVol;
    delete [] VpVol_fr;
    */
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
