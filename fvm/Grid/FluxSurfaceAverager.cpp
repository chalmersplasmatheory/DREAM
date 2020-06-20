/**
 * Implementation of the FluxSurfaceAverager class which 
 * handles everything needed to carry out flux surface 
 * averages in DREAM. Is initialized by a RadialGridGenerator
 * via the method SetReferenceMagneticFieldData.
 */

#include "FVM/Grid/FluxSurfaceAverager.hpp"

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

}

/**
 * Destructor
 */
FluxSurfaceAverager::~FluxSurfaceAverager(){
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
    real_t Bmin, Bmax, VpVol;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL){
        Bmin = rGrid->GetBmin_f(ir);
        Bmax = rGrid->GetBmax_f(ir);
        VpVol = rGrid->GetVpVol_f(ir);
    } else {
        Bmin = rGrid->GetBmin(ir);
        Bmax = rGrid->GetBmax(ir);
        VpVol = rGrid->GetVpVol(ir);
    }
    if (Bmin == Bmax){
        return F(1.0,1.0,1.0);
    } else 
        return EvaluateFluxSurfaceIntegral(ir,fluxGridType, F) / VpVol;
}


/**
 * The function is used in the evaluation of the effective passing fraction,
 * and represents x / <1-x B/Bmax>
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
    real_t Bmin;
    if (fluxGridType==FLUXGRIDTYPE_RADIAL){
        Bmin     = rGrid->GetBmin_f(ir);
    } else {
        Bmin     = rGrid->GetBmin(ir);
    }

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
    gsl_integration_fixed_free(gsl_w);
    delete [] theta;
    delete [] weights;
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
    
    std::function<real_t(real_t,real_t,real_t)>  QuadWeightFunction;
    
    const gsl_integration_fixed_type *quadratureRule;
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
            gsl_adaptive = gsl_integration_workspace_alloc(1000);
            integrateAdaptive = true;
            return;
        default:
            throw FVMException("Quadrature rule '%d' not supported by FluxSurfaceAverager.", q_method);            
    }

    if(geometryIsSymmetric)
        theta_max = M_PI;
    else 
        theta_max = 2*M_PI;

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
    real_t Bmin;
    if (fluxGridType==FLUXGRIDTYPE_RADIAL)
        Bmin = rGrid->GetBmin_f(ir);
    else
        Bmin = rGrid->GetBmin(ir);
    return evaluateXiAtB(xi0, B->evaluateAtTheta(ir,theta,fluxGridType) / Bmin ); 
}
        

