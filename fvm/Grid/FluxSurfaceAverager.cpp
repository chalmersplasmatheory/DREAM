/**
 * Implementation of the FluxSurfaceAverager class which 
 * handles everything needed to carry out flux surface 
 * averages in DREAM. 
 * 
 * It is initialized (or refreshed) with Rebuild(), which should
 * be called after RadialGridGenerator has called 
 * SetReferenceMagneticFieldData(...).
 */

#include "FVM/Grid/FluxSurfaceAverager.hpp"
#include "gsl/gsl_errno.h"
using namespace std;
using namespace DREAM::FVM;


/**
 * Constructor.
 */
FluxSurfaceAverager::FluxSurfaceAverager(
    RadialGrid *g, RadialGridGenerator *rgg, bool geometryIsSymmetric, len_t ntheta_interp,
    interp_method i_method, quadrature_method q_method
) : rGrid(g), gridGenerator(rgg), geometryIsSymmetric(geometryIsSymmetric), 
    ntheta_interp(ntheta_interp) {
    const gsl_interp_type *interpolationMethod;
    switch(i_method){
        case INTERP_LINEAR:
            interpolationMethod = gsl_interp_linear;
            break;
        case INTERP_STEFFEN:{
            if(ntheta_interp>2)
                interpolationMethod = gsl_interp_steffen;
            else
                interpolationMethod = gsl_interp_linear;
            break;
        }
        default:
            throw FVMException("Interpolation method '%d' not supported by FluxSurfaceAverager.", i_method);
    }

    InitializeQuadrature(q_method);

    BOverBmin = new FluxSurfaceQuantity(rGrid, [rgg,g](len_t ir, real_t theta){if(!g->GetBmin(ir)) return 1.0; else return rgg->BAtTheta(ir,theta)/g->GetBmin(ir);}, [rgg,g](len_t ir, real_t theta){if(!g->GetBmin_f(ir)) return 1.0; else return rgg->BAtTheta_f(ir,theta)/g->GetBmin_f(ir);}, interpolationMethod);
    Jacobian  = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta){return rgg->JacobianAtTheta(ir,theta);},           [rgg](len_t ir, real_t theta){return rgg->JacobianAtTheta_f(ir,theta);}, interpolationMethod);
    ROverR0   = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta){return rgg->ROverR0AtTheta(ir,theta);},            [rgg](len_t ir, real_t theta){return rgg->ROverR0AtTheta_f(ir,theta);}, interpolationMethod);
    NablaR2   = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta){return rgg->NablaR2AtTheta(ir,theta);},            [rgg](len_t ir, real_t theta){return rgg->NablaR2AtTheta_f(ir,theta);}, interpolationMethod);

    // Use the Brent algorithm for root finding in determining the theta bounce points
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_fsolver = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    qaws_table = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0, 0);
}

/**
 * Destructor
 */
FluxSurfaceAverager::~FluxSurfaceAverager(){
    gsl_root_fsolver_free(gsl_fsolver);

    DeallocateQuadrature();
    DeallocateReferenceData();
    delete BOverBmin;
    delete Jacobian;
    delete ROverR0;
    delete NablaR2;
}


/**
 * (Re-)Initializes everyting required to perform flux surface averages.
 * Should be called after SetReferenceMagneticFieldData(...).
 */
void FluxSurfaceAverager::Rebuild(){
    this->nr = rGrid->GetNr();

    // if using fixed quadrature, store all quantities on the theta grid
    if(!integrateAdaptive){
        BOverBmin->InterpolateMagneticDataToTheta(theta, ntheta_interp);
        Jacobian->InterpolateMagneticDataToTheta(theta, ntheta_interp);
        ROverR0->InterpolateMagneticDataToTheta(theta, ntheta_interp);
        NablaR2->InterpolateMagneticDataToTheta(theta, ntheta_interp);
    }

    // Calculate flux-surface averaged Jacobian and hand over to RadialGrid.
    function<real_t(real_t,real_t,real_t)> unityFunc 
               = [](real_t,real_t,real_t){return 1;};
    int_t unityList[4] = {0,0,0,1};
    real_t *VpVol   = new real_t[nr];
    real_t *VpVol_f = new real_t[nr+1];    
    for(len_t ir=0; ir<nr;  ir++)
        VpVol[ir]   = EvaluateFluxSurfaceIntegral(ir, FLUXGRIDTYPE_DISTRIBUTION, unityFunc, unityList);
    for(len_t ir=0; ir<=nr; ir++)
        VpVol_f[ir] = EvaluateFluxSurfaceIntegral(ir, FLUXGRIDTYPE_RADIAL, unityFunc, unityList);
    
    rGrid->SetVpVol(VpVol,VpVol_f);
}


/**
 *  Evaluates the flux surface average <F> of a function F = F(B/Bmin, R/R0, |nabla r|^2) on radial grid point ir. 
 */
real_t FluxSurfaceAverager::CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t)> F, int_t *F_list){
    real_t VpVol = GetVpVol(ir,fluxGridType);

    // treat singular point r=0 separately where orbit parameters are constant 
    if(VpVol == 0) 
        return F(1,1,1);

    // otherwise use regular method
    return EvaluateFluxSurfaceIntegral(ir,fluxGridType, F, F_list) / VpVol;
}


/**
 * The function returns the integrand of the flux surface integral at theta,
 * and is used with the adaptive quadrature
 */
struct FluxSurfaceIntegralParams {
    std::function<real_t(real_t,real_t,real_t)> Function; len_t ir; real_t Bmin;
    FluxSurfaceAverager *FSA; fluxGridType fgType;
};
real_t FluxSurfaceAverager::FluxSurfaceIntegralFunction(real_t theta, void *p){
    struct FluxSurfaceIntegralParams *params = (struct FluxSurfaceIntegralParams *) p;
    len_t ir = params->ir;
    fluxGridType fluxGridType = params->fgType;
    real_t Bmin = params->Bmin;

    real_t B,Jacobian,ROverR0,NablaR2;
    params->FSA->GeometricQuantitiesAtTheta(ir,theta,B,Jacobian,ROverR0,NablaR2,fluxGridType);
    real_t BOverBmin=1;
    if(Bmin != 0)
        BOverBmin = B/Bmin;
    std::function<real_t(real_t,real_t,real_t)> F = params->Function;    
    return 2*M_PI*Jacobian*F(BOverBmin, ROverR0, NablaR2);
}


/**
 * Core function of this class: evaluates the flux surface integral 
 *    FluxSurfaceIntegral(X) = \int J*X dphi dtheta
 * taken over toroidal and poloidal angle, weighted by the 
 * spatial jacobian J. See doc/notes/theory section on 
 * Flux surface averages for further details.
 */
real_t FluxSurfaceAverager::EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t)> F, int_t */*F_list*/){
    real_t fluxSurfaceIntegral = 0;

    // Calculate using fixed quadrature:
    if(!integrateAdaptive){
        const real_t *BOverBmin = this->BOverBmin->GetData(ir, fluxGridType);
        const real_t *Jacobian  = this->Jacobian->GetData(ir,fluxGridType);
        const real_t *ROverR0   = this->ROverR0->GetData(ir, fluxGridType);
        const real_t *NablaR2   = this->NablaR2->GetData(ir, fluxGridType);
    
        for (len_t it = 0; it<ntheta_interp; it++)
            fluxSurfaceIntegral += weights[it] * Jacobian[it] 
                * F(BOverBmin[it], ROverR0[it], NablaR2[it]);
        fluxSurfaceIntegral *= 2*M_PI;
    // or by using adaptive quadrature:
    } else {
        gsl_function GSL_func; 
        FluxSurfaceIntegralParams params = {F, ir, GetBmin(ir,fluxGridType), this, fluxGridType}; 
        GSL_func.function = &(FluxSurfaceIntegralFunction);
        GSL_func.params = &params;
        real_t epsabs = 0, epsrel = 1e-4, lim = gsl_adaptive->limit, error;
        gsl_integration_qag(&GSL_func, 0, theta_max,epsabs,epsrel,lim,QAG_KEY,gsl_adaptive,&fluxSurfaceIntegral, &error);
    }
    return fluxSurfaceIntegral;  
} 



/**
 * Deallocate quadrature.
 */
void FluxSurfaceAverager::DeallocateQuadrature(){
    gsl_integration_workspace_free(gsl_adaptive);
    gsl_integration_workspace_free(gsl_adaptive_outer);
    gsl_integration_qaws_table_free(qaws_table);
    if(gsl_w != nullptr)
        gsl_integration_fixed_free(gsl_w);
}

/**
 * Sets quantities related to the quadrature rule used for 
 * evaluating the flux surface integrals. We use a quadrature
 *      integral( w(x, xmin, xmax)*F(x), xmin, xmax )
 *          = sum_i w_i F(x_i)
 * where theta is x_i, weights is w_i and quadratureWeightFunction 
 * is w(x, xmin, xmax).  ntheta_interp is the number of points i 
 * in the sum. The same quadrature is used on all radii.
 * Refer to GSL integration documentation for possible quadrature rules 
 * and corresponding weight functions:
 *      https://www.gnu.org/software/gsl/doc/html/integration.html
 */
void FluxSurfaceAverager::InitializeQuadrature(quadrature_method q_method){
    gsl_adaptive = gsl_integration_workspace_alloc(1000);
    gsl_adaptive_outer = gsl_integration_workspace_alloc(1000);
    std::function<real_t(real_t,real_t,real_t)>  QuadWeightFunction;
    if(geometryIsSymmetric)
        theta_max = M_PI;
    else 
        theta_max = 2*M_PI;

    const gsl_integration_fixed_type *quadratureRule = nullptr;
    switch(q_method){
        case QUAD_FIXED_LEGENDRE:
            // Legendre quadrature is often best for smooth finite functions on a finite interval
            quadratureRule = gsl_integration_fixed_legendre;
            QuadWeightFunction = [](real_t /*x*/, real_t /*x_min*/, real_t /*x_max*/)
                                    {return 1;};
            break;
        case QUAD_FIXED_CHEBYSHEV:
            // Chebyshev quadrature may be better for integration along trapped orbits
            // where the metric takes the form of this weight function.
            quadratureRule = gsl_integration_fixed_chebyshev;
            QuadWeightFunction = [](real_t x, real_t x_min, real_t x_max)
                                    {return 1/sqrt((x_max-x)*(x-x_min) );};
            break;
        case QUAD_ADAPTIVE:
            integrateAdaptive = true;
            return;
        default:
            throw FVMException("Quadrature rule '%d' not supported by FluxSurfaceAverager.", q_method);            
    }

    
    gsl_w = gsl_integration_fixed_alloc(quadratureRule,ntheta_interp,0,theta_max,0,0);
    this->theta = gsl_w->x;
    this->weights = gsl_w->weights;

    // We divide the weights by the QuadWeightFunction to cast the
    // equation on the form required by the quadrature rule.
    for(len_t it=0; it<ntheta_interp; it++)
        weights[it]/= QuadWeightFunction(theta[it],0,theta_max);

    // If symmetric field we integrate from 0 to pi and multiply result by 2. 
    if(geometryIsSymmetric)
        for(len_t it=0; it<ntheta_interp; it++)
            weights[it] *= 2;
}


/**
 * Sets reference magnetic field quantities that have been
 * generated by a RadialGridGenerator (FluxSurfaceAverager
 * takes ownership of these and will deallocate them).
 */
void FluxSurfaceAverager::SetReferenceMagneticFieldData(
    real_t *theta_Bmin, real_t *theta_Bmin_f,
    real_t *theta_Bmax, real_t *theta_Bmax_f
){
    InitializeReferenceData(theta_Bmin, theta_Bmin_f, theta_Bmax, theta_Bmax_f);    
}

void FluxSurfaceAverager::InitializeReferenceData(
    real_t *theta_Bmin, real_t *theta_Bmin_f,
    real_t *theta_Bmax, real_t *theta_Bmax_f
){
    DeallocateReferenceData();
    this->theta_Bmin = theta_Bmin;
    this->theta_Bmin_f = theta_Bmin_f;
    this->theta_Bmax = theta_Bmax;
    this->theta_Bmax_f = theta_Bmax_f;    
}

void FluxSurfaceAverager::DeallocateReferenceData(){
    if(theta_Bmin==nullptr)
        return;

    delete [] theta_Bmin;
    delete [] theta_Bmin_f;
    delete [] theta_Bmax;
    delete [] theta_Bmax_f;
    
}


/**
 * Helper function to get Bmin from RadialGrid.
 */
real_t FluxSurfaceAverager::GetBmin(len_t ir,fluxGridType fluxGridType, real_t *theta_Bmin){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        if(theta_Bmin != nullptr)
            *theta_Bmin = this->theta_Bmin_f[ir];
        return rGrid->GetBmin_f(ir);
    } else {
        if(theta_Bmin != nullptr)
            *theta_Bmin = this->theta_Bmin[ir];
        return rGrid->GetBmin(ir);
    }
}
/**
 * Helper function to get Bmax from RadialGrid.
 */
real_t FluxSurfaceAverager::GetBmax(len_t ir,fluxGridType fluxGridType, real_t *theta_Bmax){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        if(theta_Bmax != nullptr)
            *theta_Bmax = this->theta_Bmax_f[ir];
        return rGrid->GetBmax_f(ir);
    } else {
        if(theta_Bmax != nullptr)
            *theta_Bmax = this->theta_Bmax[ir];
        return rGrid->GetBmax(ir);
    }
}
/**
 * Helper function to get VpVol from RadialGrid.
 */
real_t FluxSurfaceAverager::GetVpVol(len_t ir,fluxGridType fluxGridType){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL)
        return rGrid->GetVpVol_f(ir);
    else
        return rGrid->GetVpVol(ir);
}


/********************************************************************
 * The below functions evaluate bounce averages in p-xi coordinates *
 * at arbitrary (p, xi0) independently of MomentumGrids.            *
 ********************************************************************/

/** 
 * Calculations related to the bounce average of the delta function
 * appearing in the Rosenbluth-Putvinski avalanche source
 */
#include "FluxSurfaceAverager.avalancheDelta.cpp"


/**
 * Evalutes the singular product sqrt(g) * sqrt(theta-theta_b1)*sqrt(theta_b2-theta)
 * needed for the QAWS quadrature in the limits theta->theta_b1 and theta->theta_b2.
 * Is based on a Taylor expansion of 
 *      xi = sqrt(1 - B/Bmin (1-xi0^2)) 
 * at the bounce point theta = theta_b.
 */
real_t EvaluateBounceIntegrandWithQAWSWeightAtSingularity(
    len_t ir, real_t xi0, real_t theta, real_t theta_b1, real_t theta_b2, 
    real_t B, real_t Bmin, FluxSurfaceAverager *FSA, fluxGridType fgType
){
    real_t otherSqrtFactor, theta0;
    if(theta_b2-theta < theta-theta_b1){
        otherSqrtFactor = sqrt(theta - theta_b1);
        theta0 = theta_b2;
    }
    else {
        otherSqrtFactor = sqrt(theta_b2 - theta);
        theta0 = theta_b1;
    }
    real_t eps=1e-4*theta0; 
    
    real_t BPrimeAtThetaB = fabs( (FSA->BAtTheta(ir,theta0+eps,fgType) - FSA->BAtTheta(ir,theta0-eps,fgType))/(2*eps));
    real_t divergentSqrtFactorOverXiOverXi0 = sqrt(xi0*xi0/ ( (BPrimeAtThetaB/Bmin) * (1-xi0*xi0) ) );
    real_t sqrtGTimesXiOverXi0 = 2*M_PI*B/Bmin; // explicit expression for PXiMetric*xiOverXi0
    return sqrtGTimesXiOverXi0*divergentSqrtFactorOverXiOverXi0*otherSqrtFactor;
}

/**
 * Evaluates the bounce-integral integrand normalized to p^2*R0 at an arbitrary 
 * poloidal angle theta, in a form that can be used by gsl quadratures.
 */
real_t FluxSurfaceAverager::BounceIntegralFunction(real_t theta, void *par){
    struct BounceIntegralParams *params = (struct BounceIntegralParams *) par;
    
    real_t xi0 = params->xi0;
    real_t Bmin = params->Bmin;
    FluxSurfaceAverager *fluxAvg = params->fsAvg; 
    int_t *Flist_eff = params->Flist_eff;
    
    real_t B,Jacobian,ROverR0,NablaR2;
    fluxAvg->GeometricQuantitiesAtTheta(params->ir,theta,B,Jacobian,ROverR0,NablaR2,params->fgType);
    real_t BOverBmin=1;
    if(Bmin != 0)
        BOverBmin = B/Bmin;    
    real_t xiOverXi0 = MomentumGrid::evaluateXiOverXi0(xi0, BOverBmin);
    real_t Function;
    if(Flist_eff != nullptr) // if the function exponent list is provided, prioritize to build function this way
        Function = fluxAvg->AssembleBAFunc(xiOverXi0, BOverBmin, ROverR0, NablaR2,Flist_eff);
    else
        Function = params->F_eff(xiOverXi0,BOverBmin,ROverR0,NablaR2);

    // handle the singular point theta=theta_b1 or theta=theta_b2 (used by QAWS)
    if(xiOverXi0<1e-4 && params->integrateQAWS) { 
        if(!params->integrateQAWS)
            throw FVMException("FluxSurfaceAverager: Cannot evaluate the BounceInteGralFunction"
                " at theta=theta_b unless using the QAWS quadrature.");
        real_t rescaledSqrtG = EvaluateBounceIntegrandWithQAWSWeightAtSingularity(
            params->ir, xi0, theta, params->theta_b1, params->theta_b2, B, Bmin, fluxAvg, params->fgType
        );
        return 2*M_PI*Function*Jacobian*rescaledSqrtG;
    } else {
        real_t sqrtG = MomentumGrid::evaluatePXiMetricOverP2(xiOverXi0,BOverBmin);
        real_t S = 2*M_PI*Jacobian*sqrtG*Function;

        if(params->integrateQAWS) // divide by weight function in QAWS quadrature
            return S*sqrt((theta - params->theta_b1)*(params->theta_b2 - theta));
        else 
            return S;
    }
}

/**
 * Evaluates the bounce integral normalized to p^2 of the function F = F(xi/xi0, B/Bmin, ROverR0, NablaR2)  
 * at radial grid point ir and pitch xi0, using an adaptive quadrature.
 */
real_t FluxSurfaceAverager::EvaluatePXiBounceIntegralAtP(len_t ir, real_t xi0, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F, int_t *Flist){
    real_t theta_Bmin, theta_Bmax;
    real_t Bmin = GetBmin(ir,fluxGridType, &theta_Bmin);
    real_t Bmax = GetBmax(ir,fluxGridType, &theta_Bmax);
    real_t BminOverBmax;
    if(Bmin==Bmax) // handles Bmax=0 case
        BminOverBmax = 1; 
    else
        BminOverBmax = Bmin/Bmax;        

    bool isTrapped = ( (1-xi0*xi0) > BminOverBmax);
    real_t theta_b1, theta_b2;
    bool integrateQAWS = false;

    std::function<real_t(real_t,real_t,real_t,real_t)> F_eff;
    int_t Flist_copy[5];
    int_t *Flist_eff = nullptr;
    if(Flist != nullptr){
        for(len_t k=0;k<5;k++)
            Flist_copy[k] = Flist[k];
        Flist_eff = Flist_copy;
    }

    // If trapped, adds contribution from -xi0
    if (isTrapped){
        // negative-pitch particles do not exist independently; are described by the positive pitch counterpart
        if(xi0<=0)
            return 0;
        // if odd function in xi, set bounce integral to zero, otherwise multiply by 2
        if(Flist != nullptr){
            if(Flist[0]%2==1)
                return 0;
            else 
                Flist_eff[4] *= 2;
        }
        F_eff = [&](real_t x, real_t  y, real_t z, real_t w){return F(x,y,z,w) + F(-x,y,z,w);};
        FindBouncePoints(ir, Bmin, theta_Bmin, theta_Bmax, this, xi0, fluxGridType, &theta_b1, &theta_b2,gsl_fsolver,geometryIsSymmetric);
        if(theta_b1==theta_b2)
            return 0;
        if(F_eff(0,1,1,1) !=0)
            integrateQAWS = true;
    } else { 
        F_eff = F;
        theta_b1 = 0;
        theta_b2 = 2*M_PI;
    }

    gsl_function GSL_func;
    GSL_func.function = &(BounceIntegralFunction);
    BounceIntegralParams params = {
        ir,xi0,theta_b1,theta_b2,fluxGridType,
        Bmin,F_eff,Flist_eff,this,integrateQAWS
    };
    GSL_func.params = &params;
    real_t bounceIntegral, error; 

    real_t epsabs = 0, epsrel = 1e-3, lim = gsl_adaptive->limit; 
    if(integrateQAWS)
        gsl_integration_qaws(&GSL_func,theta_b1,theta_b2,qaws_table,epsabs,epsrel,lim,gsl_adaptive,&bounceIntegral,&error);
    else
        gsl_integration_qag(&GSL_func,theta_b1,theta_b2,epsabs,epsrel,lim,QAG_KEY,gsl_adaptive,&bounceIntegral,&error);
    
    return bounceIntegral;
}

// Evaluates the bounce average {F} of a function F = F(xi/xi0, B/Bmin) on grid point (ir,i,j). 
real_t FluxSurfaceAverager::CalculatePXiBounceAverageAtP(
    len_t ir, real_t xi0, fluxGridType fluxGridType, 
    std::function<real_t(real_t,real_t,real_t,real_t)> F, int_t *F_list
){    
    std::function<real_t(real_t,real_t,real_t,real_t)> FUnity = [](real_t,real_t,real_t,real_t){return 1;};
    int_t unityList[5] = {0,0,0,0,1};
    
    real_t BI = EvaluatePXiBounceIntegralAtP(ir, xi0, fluxGridType, F, F_list);
    if(!BI)
        return 0;
    real_t Vp = EvaluatePXiBounceIntegralAtP(ir, xi0, fluxGridType, FUnity, unityList);    
    /**
     * Treat special case: 
     * Either r=0 (particle doesn't move in the poloidal plane) or 
     * xi0=0 in inhomogeneous magnetic field (infinitely deeply trapped on low-field side).
     * TODO: consider more carefully what happens with xi/xi0 in the xi0=0 case.
     */
//    if(Vp == 0)  
//        return F(1,1,1,1);

    return BI / Vp;
}



/**
 * Returns the function to be bounce averaged, evaluated at poloidal angle theta.
 * If Flist is provided, returns the function
 *      BA_Func = Flist[4] * xiOverXi0^Flist[0] * BOverBmin^Flist[1] 
 *               * ROverR0^Flist[2] * NablaR2^Flist[3],
 * assuming Flist[0-3] to be integers.
 */
real_t FluxSurfaceAverager::AssembleBAFunc(real_t xiOverXi0,real_t BOverBmin, real_t ROverR0, real_t NablaR2, int_t *Flist){
    real_t BA_Func = Flist[4];
    if(Flist[0]>0)
        for(int_t k=0; k<Flist[0]; k++)
            BA_Func *= xiOverXi0;
    else if(Flist[0]<0)
        for(int_t k=0; k<-Flist[0]; k++)
            BA_Func /= xiOverXi0;
    if(Flist[1]>0)
        for(int_t k=0; k<Flist[1]; k++)
            BA_Func *= BOverBmin;
    else if(Flist[1]<0)
        for(int_t k=0; k<-Flist[1]; k++)
            BA_Func /= BOverBmin;
    if(Flist[2]>0)
        for(int_t k=0; k<Flist[2]; k++)
            BA_Func *= ROverR0;
    else if(Flist[2]<0)
        for(int_t k=0; k<-Flist[2]; k++)
            BA_Func /= ROverR0;
    if(Flist[3]>0)
        for(int_t k=0; k<Flist[3]; k++)
            BA_Func *= NablaR2;
    else if(Flist[3]<0)
        for(int_t k=0; k<-Flist[3]; k++)
            BA_Func /= NablaR2;
    return BA_Func;
}



/**
 * Returns xi(xi0)^2 and is used in the determination of the bounce points. 
 */
struct xiFuncParams {real_t xi0; len_t ir; real_t Bmin; FluxSurfaceAverager *FSA; fluxGridType fgType; };
real_t FluxSurfaceAverager::xiParticleFunction(real_t theta, void *p){
    struct xiFuncParams *params = (struct xiFuncParams *) p;    
    real_t B, Jacobian, ROverR0, NablaR2;
    params->FSA->GeometricQuantitiesAtTheta(params->ir, theta, B, Jacobian, ROverR0, NablaR2, params->fgType);
    real_t Bmin = params->Bmin;
    real_t BOverBmin = 1;
    if(Bmin)
        BOverBmin = B/Bmin;

    real_t xi0 = params->xi0; 
    return 1 - (1-xi0*xi0) * BOverBmin;
}

/**
 * Returns poloidal bounce points theta_bounce1 and theta_bounce2 with a root finding algorithm
 */
void FluxSurfaceAverager::FindBouncePoints(
    len_t ir, real_t Bmin, real_t theta_Bmin, real_t theta_Bmax, 
    FluxSurfaceAverager *FSA, real_t xi0, fluxGridType fluxGridType,  
    real_t *theta_b1, real_t *theta_b2, gsl_root_fsolver* gsl_fsolver, bool isSymmetric
){
    // define GSL function xi_particle as function of theta which evaluates B_interpolator(_fr) at theta.
    xiFuncParams xi_params = {xi0,ir,Bmin,FSA,fluxGridType}; 
    gsl_function gsl_func;
    gsl_func.function = &(xiParticleFunction);
    gsl_func.params = &xi_params;

    FluxSurfaceAverager::FindThetas(theta_Bmin, theta_Bmax, theta_b1, theta_b2, gsl_func, gsl_fsolver, isSymmetric);
}


/**
 * Finds the two poloidal angles theta1 and theta2 that solves 
 *  gsl_func.function(theta) = 0,
 * returning solutions that satisfy 
 *  gsl_function.function > 0. 
 * theta2 is assumed to exist on the interval [theta_Bmin, theta_Bmax], 
 * and theta1 on the interval [theta_Bmax-2*pi, theta_Bmin]
 */
void FluxSurfaceAverager::FindThetas(
    real_t theta_Bmin, real_t theta_Bmax, real_t *theta1, real_t *theta2, 
    gsl_function gsl_func, gsl_root_fsolver *gsl_fsolver, bool isSymmetric
){
    real_t root=0;

    real_t x_lower = theta_Bmin;
    real_t x_upper = theta_Bmax;
    FindRoot(&x_lower, &x_upper, &root, gsl_func, gsl_fsolver);

    if( gsl_func.function(x_lower, gsl_func.params) >= 0 )
        *theta2 = x_lower;
    else if ( gsl_func.function(x_upper, gsl_func.params) >= 0 )
        *theta2 = x_upper;
    else
        throw FVMException("FluxSurfaceAverager: unable to find valid theta root.");

    if(isSymmetric){
        // handle special case where global minimum is for theta>0 in symmetric fields;
        // then there can be tricky solutions in the interval [0,theta_Bmin].
        if(theta_Bmin>0 && gsl_func.function(0, gsl_func.params) < 0){
            // look for the solution in the interval [0, theta_Bmin]
            x_lower = 0;
        // otherwise the solution is mirrored
        } else {
            *theta1 = -*theta2;
            return;
        }
    } else // if not symmetric, look for solution in the remaining interval
        x_lower = theta_Bmax-2*M_PI;

    x_upper = theta_Bmin;
    FindRoot(&x_lower, &x_upper, &root, gsl_func, gsl_fsolver);

    if( gsl_func.function(x_lower, gsl_func.params) >= 0 )
        *theta1 = x_lower;
    else if( gsl_func.function(x_upper, gsl_func.params) >= 0 )
        *theta1 = x_upper;
    else
        throw FVMException("FluxSurfaceAverager: unable to find valid theta root.");    
}


/**
 * Finds the theta that solves gsl_func = 0.
 * Takes a theta interval theta \in [x_lower, x_upper] and iterates at most 
 * max_iter times (or to a relative error of rel_error) to find an estimate 
 * for theta \in [x_lower, x_upper] 
 */
void FluxSurfaceAverager::FindRoot(
    real_t *x_lower, real_t *x_upper, real_t *root, 
    gsl_function gsl_func, gsl_root_fsolver *gsl_fsolver
){
    gsl_root_fsolver_set (gsl_fsolver, &gsl_func, *x_lower, *x_upper); // finds root in [0,pi] using GSL_rootsolver_type algorithm
    int status;
    real_t epsrel = 1e-5;
    real_t epsabs = 1e-5;
    len_t max_iter = 50;    
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status   = gsl_root_fsolver_iterate (gsl_fsolver);
        *root    = gsl_root_fsolver_root (gsl_fsolver);
        *x_lower = gsl_root_fsolver_x_lower (gsl_fsolver);
        *x_upper = gsl_root_fsolver_x_upper (gsl_fsolver);
        status   = gsl_root_test_interval (*x_lower, *x_upper,
                                            epsabs, epsrel);

        if (status == GSL_SUCCESS)
            break;
    }
}


///////////////////////////////////////////////////////////////
// METHODS FOR AVERAGING A BOUNCE INTEGRAL OVER A PITCH CELL //
///////////////////////////////////////////////////////////////

const real_t NEAR_TRAPPED_BOUNDARY_THRESHOLD = 0.01; // distance in xi from trapped boundary within which we should carefully integrate fluxes
/**
 * Returns `true` if a cell with pitch faces xi_lower and xi_upper satisfies certain constraints
 * that demands a careful cell average in order to retain precision in the Vp calculation 
 * (most notably on the trapped-passing boundary where Vp has a logarithmic singularity).  
 */
bool FluxSurfaceAverager::shouldCellAverageBounceIntegral(len_t ir, real_t xi_lower, real_t xi_upper, fluxGridType fluxGridType){
    real_t xiT;
    if(fluxGridType==FLUXGRIDTYPE_RADIAL)
        xiT = rGrid->GetXi0TrappedBoundary_fr(ir);
    else
        xiT = rGrid->GetXi0TrappedBoundary(ir);
    // if the cell contains, or is near, the trapped-passing boundaries, 
    // perform more careful (and computationally expensive) cell average
    real_t d1=fabs(xi_lower)-xiT;
    real_t d2=fabs(xi_upper)-xiT;
    bool containsTrappedBoundary = (xi_lower<xiT && xi_upper >= xiT) || (xi_lower<=-xiT && xi_upper >-xiT);
    bool isNearTrappedBoundary = min(fabs(d1),fabs(d2)) < NEAR_TRAPPED_BOUNDARY_THRESHOLD;
    return containsTrappedBoundary || isNearTrappedBoundary;
}


/** 
 * Helper function for gsl integration: returns bounce integral function at xi
 */
struct PXiIntegralParams {len_t ir; fluxGridType fgType; function<real_t(real_t,real_t,real_t,real_t)> F; int_t *Flist; FluxSurfaceAverager *FSA;};
real_t FluxSurfaceAverager::evaluatePXiBounceIntegralAtXi(real_t xi0, void *par){
    PXiIntegralParams *params = (struct PXiIntegralParams*)par;
    return params->FSA->EvaluatePXiBounceIntegralAtP(params->ir, xi0, params->fgType, params->F, params->Flist);
}


/**
 * Averages the bounce integral over xi from xi_l to xi_u (which are assumed to satisfy xi_l<=xi_u)
 */
real_t FluxSurfaceAverager::EvaluateCellAveragedBounceIntegralOverP2(len_t ir, real_t xi_l, real_t xi_u, fluxGridType fluxGridType, function<real_t(real_t,real_t,real_t,real_t)> F, int_t *F_list){
    if( fabs(xi_u-xi_l) < 100*std::numeric_limits<real_t>::epsilon()) // simply evaluate the bounce integral at the point xi_u=xi_l
        return EvaluatePXiBounceIntegralAtP(ir,xi_u,fluxGridType,F,F_list);
    real_t dxi = xi_u - xi_l;
    real_t 
        partResult1 = 0,
        partResult2 = 0,
        partResult3 = 0;
    
    int key = GSL_INTEG_GAUSS41;
    real_t 
        epsabs = 0,
        epsrel = 1e-3,
        lim = gsl_adaptive->limit,
        error;
    
    real_t xiT;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        xiT = rGrid->GetXi0TrappedBoundary_fr(ir);
    else 
        xiT = rGrid->GetXi0TrappedBoundary(ir);
    PXiIntegralParams params = {ir, fluxGridType, F, F_list, this}; 
    gsl_function gsl_func;
    gsl_func.function = &(evaluatePXiBounceIntegralAtXi);
    gsl_func.params = &params;
    
    // contribution from negative pitch passing region
    if(xi_l < -xiT){
        real_t pts[2] = {xi_l,-xiT};
        int npts = 2;
        if(xi_u < -xiT)
            gsl_integration_qag(&gsl_func,xi_l,xi_u,epsabs,epsrel,lim,key,gsl_adaptive_outer,&partResult1,&error);
        else{}
//            gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive_outer,&partResult1,&error);
    }
    // contribution from positive pitch trapped region
    if(xi_u>0 && xi_l<xiT){
        if(xi_u<xiT){
            if(xi_l<0)
                gsl_integration_qag(&gsl_func,0,xi_u,epsabs,epsrel,lim,key,gsl_adaptive_outer,&partResult2,&error);
            else 
                gsl_integration_qag(&gsl_func,xi_l,xi_u,epsabs,epsrel,lim,key,gsl_adaptive_outer,&partResult2,&error);
        } else {
            if(xi_l<0){
                real_t pts[2] = {0,xiT};
                int npts = 2;
                gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive_outer,&partResult2,&error);
            } else {
                real_t pts[2] = {xi_l,xiT};
                int npts = 2;
                gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive_outer,&partResult2,&error);
            }
        }
    }
    // contribution from positive pitch passing region
    if(xi_u>xiT){
        real_t pts[2] = {xiT,xi_u};
        int npts = 2;
        if(xi_l<=xiT)
            gsl_integration_qagp(&gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive_outer,&partResult3,&error);
        else 
            gsl_integration_qag(&gsl_func,xi_l,xi_u,epsabs,epsrel,lim,key,gsl_adaptive_outer,&partResult3,&error);
    }    

    return (partResult1+partResult2+partResult3)/dxi;
}
