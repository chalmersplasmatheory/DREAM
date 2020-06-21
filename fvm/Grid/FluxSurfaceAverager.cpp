/**
 * Implementation of the FluxSurfaceAverager class which 
 * handles everything needed to carry out flux surface 
 * averages in DREAM. Is initialized by a RadialGridGenerator
 * via the method SetReferenceMagneticFieldData.
 */

#include "FVM/Grid/FluxSurfaceAverager.hpp"
#include "gsl/gsl_errno.h"
using namespace std;
using namespace DREAM::FVM;


/**
 * Constructor.
 */
FluxSurfaceAverager::FluxSurfaceAverager(
    RadialGrid *g, bool geometryIsSymmetric, len_t ntheta_interp,
    interp_method i_method, quadrature_method q_method
) : rGrid(g), geometryIsSymmetric(geometryIsSymmetric), 
    ntheta_interp(ntheta_interp) {
    const gsl_interp_type *interpolationMethod;
    switch(i_method){
        case INTERP_LINEAR:
            interpolationMethod = gsl_interp_linear;
            break;
        case INTERP_STEFFEN:
            interpolationMethod = gsl_interp_steffen;
            break;
        default:
            throw FVMException("Interpolation method '%d' not supported by FluxSurfaceAverager.", i_method);
    }

    InitializeQuadrature(q_method);

    B        = new FluxSurfaceQuantity(rGrid, interpolationMethod);
    Jacobian = new FluxSurfaceQuantity(rGrid, interpolationMethod);
    ROverR0  = new FluxSurfaceQuantity(rGrid, interpolationMethod);
    NablaR2  = new FluxSurfaceQuantity(rGrid, interpolationMethod);

    // Use the Brent algorithm for root finding in determining the theta bounce points
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_fsolver = gsl_root_fsolver_alloc (GSL_rootsolver_type);

}

/**
 * Destructor
 */
FluxSurfaceAverager::~FluxSurfaceAverager(){
    gsl_root_fsolver_free(gsl_fsolver);

    DeallocateQuadrature();
    delete B;
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


    if(!integrateAdaptive){
        B->InterpolateMagneticDataToTheta(theta, ntheta_interp);
        Jacobian->InterpolateMagneticDataToTheta(theta, ntheta_interp);
        ROverR0->InterpolateMagneticDataToTheta(theta, ntheta_interp);
        NablaR2->InterpolateMagneticDataToTheta(theta, ntheta_interp);
    }
    function<real_t(real_t,real_t,real_t)> unityFunc 
               = [](real_t,real_t,real_t){return 1;};
    real_t *VpVol   = new real_t[nr];
    real_t *VpVol_f = new real_t[nr+1];    
    for(len_t ir=0; ir<nr;  ir++)
        VpVol[ir]   = EvaluateFluxSurfaceIntegral(ir, FLUXGRIDTYPE_DISTRIBUTION, unityFunc);
    for(len_t ir=0; ir<=nr; ir++)
        VpVol_f[ir] = EvaluateFluxSurfaceIntegral(ir, FLUXGRIDTYPE_RADIAL, unityFunc);
    
    rGrid->SetVpVol(VpVol,VpVol_f);
}


// Evaluates the flux surface average <F> of a function F = F(B/Bmin, R/R0, |nabla r|^2) on radial grid point ir. 
real_t FluxSurfaceAverager::CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t)> F){
    real_t Bmin = GetBmin(ir,fluxGridType);
    real_t VpVol = GetVpVol(ir,fluxGridType);
    real_t Bmax;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        Bmax = rGrid->GetBmax_f(ir);
    else
        Bmax = rGrid->GetBmax(ir);
    if (Bmin == Bmax){
        return F(1.0,1.0,1.0);
    } else 
        return EvaluateFluxSurfaceIntegral(ir,fluxGridType, F) / VpVol;
}


/**
 * The function returns the integrand of the flux surface integral at theta
 */
struct FluxSurfaceIntegralParams {
    std::function<real_t(real_t,real_t,real_t)> Function; len_t ir; real_t Bmin;
    FluxSurfaceQuantity *B; FluxSurfaceQuantity *Jacobian; FluxSurfaceQuantity *ROverR0; 
    FluxSurfaceQuantity *NablaR2; fluxGridType fgType;
};
real_t FluxSurfaceAverager::FluxSurfaceIntegralFunction(real_t theta, void *p){
    struct FluxSurfaceIntegralParams *params = (struct FluxSurfaceIntegralParams *) p;
    len_t ir = params->ir;
    fluxGridType fluxGridType = params->fgType;
    real_t B = params->B->evaluateAtTheta(ir, theta, fluxGridType);
    real_t Jacobian = params->Jacobian->evaluateAtTheta(ir, theta, fluxGridType);
    real_t ROverR0 = params->ROverR0->evaluateAtTheta(ir, theta, fluxGridType);
    real_t NablaR2 = params->NablaR2->evaluateAtTheta(ir, theta, fluxGridType);
    real_t Bmin = params->Bmin;
    std::function<real_t(real_t,real_t,real_t)> F = params->Function;

    return 2*M_PI*Jacobian*F(B/Bmin, ROverR0, NablaR2);
}


real_t FluxSurfaceAverager::EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t)> F){
    real_t fluxSurfaceIntegral = 0;
    real_t Bmin = GetBmin(ir,fluxGridType);

    if(!integrateAdaptive){
        const real_t *B = this->B->GetData(ir, fluxGridType);
        const real_t *Jacobian = this->Jacobian->GetData(ir,fluxGridType);
        const real_t *ROverR0 = this->ROverR0->GetData(ir, fluxGridType);
        const real_t *NablaR2 = this->NablaR2->GetData(ir, fluxGridType);
        
    
        for (len_t it = 0; it<ntheta_interp; it++){
            fluxSurfaceIntegral += 2*M_PI*weights[it] * Jacobian[it] 
                * F(B[it]/Bmin, ROverR0[it], NablaR2[it]);
        }
    } else {
        gsl_function GSL_func; 
        FluxSurfaceIntegralParams params = {F, ir, Bmin, B, Jacobian, ROverR0, NablaR2, fluxGridType}; 
        GSL_func.function = &(FluxSurfaceIntegralFunction);
        GSL_func.params = &params;
        real_t epsabs = 0, epsrel = 1e-4, lim = gsl_adaptive->limit, error;
//        len_t key = GSL_INTEG_GAUSS31; // intermediate-order method 
        gsl_integration_qags(&GSL_func, 0, theta_max,epsabs,epsrel,lim,gsl_adaptive,&fluxSurfaceIntegral, &error);
    }
    return fluxSurfaceIntegral;  
} 



/**
 * Deallocate quadrature.
 */
void FluxSurfaceAverager::DeallocateQuadrature(){
    gsl_integration_workspace_free(gsl_adaptive);

    if(gsl_w != nullptr){
        gsl_integration_fixed_free(gsl_w);
    }
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
    std::function<real_t(real_t,real_t,real_t)>  QuadWeightFunction;
    if(geometryIsSymmetric)
        theta_max = M_PI;
    else 
        theta_max = 2*M_PI;

    const gsl_integration_fixed_type *quadratureRule = nullptr;
    switch(q_method){
        case QUAD_FIXED_LEGENDRE:
            // Legendre quadrature is often good for smooth finite functions on a finite interval
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
    len_t ntheta_ref, real_t *theta_ref,
    real_t **B_ref, real_t **B_ref_f,
    real_t **Jacobian_ref, real_t **Jacobian_ref_f,
    real_t **ROverR0_ref ,real_t **ROverR0_ref_f, 
    real_t **NablaR2_ref, real_t **NablaR2_ref_f
){
    B->Initialize(B_ref,B_ref_f, theta_ref, ntheta_ref);
    Jacobian->Initialize(Jacobian_ref, Jacobian_ref_f, theta_ref, ntheta_ref);
    ROverR0->Initialize( ROverR0_ref,  ROverR0_ref_f,  theta_ref, ntheta_ref);
    NablaR2->Initialize( NablaR2_ref,  NablaR2_ref_f,  theta_ref, ntheta_ref);
}


 real_t FluxSurfaceAverager::evaluateXiAtTheta(len_t ir, real_t xi0, real_t theta, fluxGridType fluxGridType){ 
    return evaluateXiAtB(xi0, B->evaluateAtTheta(ir,theta,fluxGridType) / GetBmin(ir,fluxGridType) ); 
}


real_t FluxSurfaceAverager::GetBmin(len_t ir,fluxGridType fluxGridType){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        return rGrid->GetBmin_f(ir);
    } else {
        return rGrid->GetBmin(ir);
    }
}
real_t FluxSurfaceAverager::GetVpVol(len_t ir,fluxGridType fluxGridType){
    if (fluxGridType == FLUXGRIDTYPE_RADIAL){
        return rGrid->GetVpVol_f(ir);
    } else {
        return rGrid->GetVpVol(ir);
    }
}


/**
 * Remainder: calculation of bounce average at arbitrary (p,xi) 
 * (independent of momentum grid)
 */


struct generalBounceIntegralParams {
    len_t ir; real_t xi0; real_t p; fluxGridType fgType; 
    real_t Bmin; std::function<real_t(real_t,real_t,real_t,real_t)> F_eff; 
    FluxSurfaceAverager *fsAvg;};
real_t generalBounceIntegralFunc(real_t theta, void *par){
    struct generalBounceIntegralParams *params = (struct generalBounceIntegralParams *) par;
    
    len_t ir = params->ir;
    real_t xi0 = params->xi0;
    real_t p = params->p;
    fluxGridType fluxGridType = params->fgType;
    real_t Bmin = params->Bmin;
    FluxSurfaceAverager *fluxAvg = params->fsAvg; 
    std::function<real_t(real_t,real_t,real_t,real_t)> F_eff = params->F_eff;
    real_t B = fluxAvg->GetB()->evaluateAtTheta(ir,theta,fluxGridType); 
    real_t Jacobian = fluxAvg->GetJacobian()->evaluateAtTheta(ir,theta,fluxGridType);
    real_t ROverR0 = fluxAvg->GetROverR0()->evaluateAtTheta(ir,theta,fluxGridType);
    real_t NablaR2 = fluxAvg->GetNablaR2()->evaluateAtTheta(ir,theta,fluxGridType);
    real_t sqrtG = MomentumGrid::evaluatePXiMetricOverP2(p,xi0,B,Bmin);
    real_t xi0Sq = xi0*xi0;
    real_t xiOverXi0 = sqrt( (1-B/Bmin*(1-xi0Sq))/xi0Sq );
    real_t F =  F_eff(xiOverXi0,B/Bmin,ROverR0,NablaR2);
    return 2*M_PI*Jacobian*sqrtG*F;
}

// Evaluates the bounce integral of the function F = F(xi/xi0, B/Bmin) at  
// radial grid point ir, momentum p and pitch xi0, using an adaptive quadrature.
real_t FluxSurfaceAverager::EvaluatePXiBounceIntegralAtP(len_t ir, real_t p, real_t xi0, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    real_t Bmin = GetBmin(ir,fluxGridType);
    real_t Bmax;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        Bmax = rGrid->GetBmax_f(ir);
    else
        Bmax = rGrid->GetBmax(ir);
    
    if(xi0*xi0 < 1e-30){
        return 0;
    }
    std::function<real_t(real_t,real_t,real_t,real_t)> F_eff;
    bool isTrapped = (Bmax/Bmin * (1-xi0*xi0) > 1);
    // If trapped, adds contribution from -xi0, since negative xi0 are presumably not kept on the grid.
    real_t theta_b1, theta_b2;
    if (isTrapped){
        F_eff = [&](real_t x, real_t  y, real_t z, real_t w){return  F(x,y,z,w) + F(-x,y,z,w) ;};
        FindBouncePoints(ir, Bmin, this->B, xi0, fluxGridType, &theta_b1, &theta_b2,gsl_fsolver, geometryIsSymmetric);
    } else { 
        F_eff = F;
        theta_b1 = 0;
        theta_b2 = 2*M_PI;
    }

    gsl_function GSL_func;
    GSL_func.function = &(generalBounceIntegralFunc);
    generalBounceIntegralParams params = {ir,xi0,p,fluxGridType,Bmin,F_eff,this};
    GSL_func.params = &params;
    real_t bounceIntegral, error; 

    real_t epsabs = 0, epsrel = 1e-4, lim = gsl_adaptive->limit; 
    gsl_integration_qags(&GSL_func,theta_b1,theta_b2,epsabs,epsrel,lim,gsl_adaptive,&bounceIntegral,&error);
    return bounceIntegral;
}


real_t singularBounceAverageFunc(real_t x, void *par){
    struct generalBounceIntegralParams *params = (struct generalBounceIntegralParams *) par;
    std::function<real_t(real_t,real_t,real_t,real_t)> F_eff = params->F_eff;
    return F_eff(sqrt(1-x*x),1,1,1)/sqrt(1-x*x);
}
// Evaluates the bounce average {F} of a function F = F(xi/xi0, B/Bmin) on grid point (ir,i,j). 
real_t FluxSurfaceAverager::CalculatePXiBounceAverageAtP(len_t ir, real_t p, real_t xi0, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F){
    
    real_t Bmin = GetBmin(ir,fluxGridType);

    if(xi0*xi0 < 1e-30) {
        generalBounceIntegralParams params = {ir,xi0,p,fluxGridType,Bmin,[F](real_t x, real_t y, real_t z, real_t w){return F(x,y,z,w)+F(-x,y,z,w);},this}; 
        real_t BounceAverageSingular, error;
        gsl_function GSL_func;
        GSL_func.function = &(singularBounceAverageFunc);
        GSL_func.params = &params;
        gsl_integration_qags(&GSL_func,0,1,0,1e-4,gsl_adaptive->limit,gsl_adaptive,&BounceAverageSingular,&error);
        return BounceAverageSingular / M_PI;
    } else {
        std::function<real_t(real_t,real_t,real_t,real_t)> FUnity = [](real_t,real_t,real_t,real_t){return 1;};
        return EvaluatePXiBounceIntegralAtP(ir, p, xi0, fluxGridType, F) / EvaluatePXiBounceIntegralAtP(ir, p, xi0, fluxGridType, FUnity);
    }
}


/**
 * The function is used in the evaluation of the effective passing fraction,
 * and represents x / <1-x B/Bmax>
 */
struct xiFuncParams {real_t xi0; len_t ir; real_t Bmin; const FluxSurfaceQuantity *B; fluxGridType fgType; };
real_t FluxSurfaceAverager::xiParticleFunction(real_t theta, void *p){
    struct xiFuncParams *params = (struct xiFuncParams *) p;
    
    real_t xi0 = params->xi0; 
    len_t ir   = params->ir;
    fluxGridType fluxGridType = params->fgType;
    real_t Bmin = params->Bmin;
    const FluxSurfaceQuantity *B = params->B;
    return 1 - (1-xi0*xi0) *B->evaluateAtTheta(ir,theta,fluxGridType) / Bmin;
}

// calculates theta_bounce1 and theta_bounce2 with a root finding algorithm
void FluxSurfaceAverager::FindBouncePoints(len_t ir, real_t Bmin, const FluxSurfaceQuantity *B, real_t xi0, fluxGridType fluxGridType,  real_t *theta_b1, real_t *theta_b2, gsl_root_fsolver* gsl_fsolver, bool geometryIsSymmetric){
    // define GSL function xi_particle as function of theta which evaluates B_interpolator(_fr) at theta.
    xiFuncParams xi_params = {xi0,ir,Bmin,B,fluxGridType}; 
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
    FindThetaBounceRoots(&x_lower, &x_upper, &root, gsl_func, gsl_fsolver);
    
    // In symmetric field, this root corresponds to theta_b2
    if (geometryIsSymmetric){
        *theta_b2 = x_lower;
        *theta_b1 = -x_lower;
    } else {

        // if xi(theta = x_upper + epsilon) is real, the root 
        // corresponds to the lower bounce point theta_b1 
        if ( xiParticleFunction(x_upper,&xi_params) > 0 ){
            *theta_b1 = x_upper;
            x_lower = M_PI; 
            x_upper = 2*M_PI;
            FindThetaBounceRoots(&x_lower, &x_upper,&root, gsl_func,gsl_fsolver);
            *theta_b2 = x_lower;
        } else {
            *theta_b2 = x_lower;
            x_lower = -M_PI;
            x_upper = 0;
            FindThetaBounceRoots(&x_lower, &x_upper, &root, gsl_func,gsl_fsolver);
            *theta_b1 = x_upper; 
        }
    }
}



// Takes a theta interval theta \in [x_lower, x_upper] and iterates at most max_iter (=15) times (or to a relative error of 0.001)
// to find an estimate for the bounce point theta_bounce \in [x_lower, x_upper]
void FluxSurfaceAverager::FindThetaBounceRoots(real_t *x_lower, real_t *x_upper, real_t *root, gsl_function gsl_func, gsl_root_fsolver *gsl_fsolver){
    gsl_root_fsolver_set (gsl_fsolver, &gsl_func, *x_lower, *x_upper); // finds root in [0,pi] using GSL_rootsolver_type algorithm
    int status;
    real_t rel_error = 1e-5;
    len_t max_iter = 15;    
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status   = gsl_root_fsolver_iterate (gsl_fsolver);
        *root    = gsl_root_fsolver_root (gsl_fsolver);
        *x_lower = gsl_root_fsolver_x_lower (gsl_fsolver);
        *x_upper = gsl_root_fsolver_x_upper (gsl_fsolver);
        status   = gsl_root_test_interval (*x_lower, *x_upper,
                                            0, rel_error);

        if (status == GSL_SUCCESS){
            break;
        }
    }
}