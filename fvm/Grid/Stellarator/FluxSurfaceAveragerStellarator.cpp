/**
 * Implementation of the FluxSurfaceAveragerStellarator class which 
 * handles everything needed to carry out flux surface 
 * averages in DREAM. 
 * 
 * It is initialized (or refreshed) with Rebuild(), which should
 * be called after RadialGridGeneratorStellarator has called 
 * SetReferenceMagneticFieldData(...).
 */

#include "FVM/Grid/Stellarator/FluxSurfaceAveragerStellarator.hpp"
#include "gsl/gsl_errno.h"
using namespace std;
using namespace DREAM::FVM;


/**
 * Constructor.
 */
FluxSurfaceAveragerStellarator::FluxSurfaceAveragerStellarator(
    RadialGridStellarator *g, RadialGridGeneratorStellarator *rgg, len_t nfp, len_t ntheta_interp, len_t nphi_interp, 
    interp_method i_method, quadrature_method q_method
) : FluxSurfaceAverager(g, rgg, false, ntheta_interp), 
    rGrid(g), gridGenerator(rgg), nfp(nfp), 
    nphi_interp(nphi_interp) {
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
            throw FVMException("Interpolation method '%d' not supported by FluxSurfaceAveragerStellarator.", i_method);
    }

    InitializeQuadrature(q_method);

    BOverBmin   = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta, real_t phi){return rgg->BAtThetaPhi(ir,theta, phi);},            [rgg](len_t ir, real_t theta, real_t phi){return rgg->BAtThetaPhi_f(ir,theta, phi);}, interpolationMethod);
    Jacobian    = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta, real_t phi){return rgg->JacobianAtThetaPhi(ir,theta, phi);},            [rgg](len_t ir, real_t theta, real_t phi){return rgg->JacobianAtThetaPhi_f(ir,theta, phi);}, interpolationMethod);
    BdotGradphi = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta, real_t phi){return rgg->BdotGradphiAtThetaPhi(ir,theta, phi);},            [rgg](len_t ir, real_t theta, real_t phi){return rgg->BdotGradphiAtThetaPhi_f(ir,theta, phi);}, interpolationMethod);
    gttOverJ2   = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta, real_t phi){return rgg->gttAtThetaPhi(ir,theta, phi);},            [rgg](len_t ir, real_t theta, real_t phi){return rgg->gttAtThetaPhi_f(ir,theta, phi);}, interpolationMethod);
    gtpOverJ2   = new FluxSurfaceQuantity(rGrid, [rgg](len_t ir, real_t theta, real_t phi){return rgg->gtpAtThetaPhi(ir,theta, phi);},            [rgg](len_t ir, real_t theta, real_t phi){return rgg->gtpAtThetaPhi_f(ir,theta, phi);}, interpolationMethod);
    // TODO: Remove when BA?
    ROverR0     = new FluxSurfaceQuantity(rGrid, [rgg](len_t, real_t){return 1.;},            [rgg](len_t, real_t){return 1.;}, interpolationMethod);
    NablaR2     = new FluxSurfaceQuantity(rGrid, [rgg](len_t, real_t){return 1.;},            [rgg](len_t, real_t){return 1.;}, interpolationMethod);

    // Use the Brent algorithm for root finding in determining the theta bounce points
    // TODO: Take back if BA
    //const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    //gsl_fsolver = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    //qaws_table = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0, 0);
}

/**
 * Destructor
 */
FluxSurfaceAveragerStellarator::~FluxSurfaceAveragerStellarator(){
    // TODO: Take back if BA
    //gsl_root_fsolver_free(gsl_fsolver);

    DeallocateQuadrature();
    DeallocateReferenceData();
    delete BOverBmin;
    delete Jacobian;
    delete BdotGradphi;
    delete gttOverJ2;
    delete gtpOverJ2;
    delete ROverR0;
    delete NablaR2;
}


/**
 * (Re-)Initializes everyting required to perform flux surface averages.
 * Should be called after SetReferenceMagneticFieldData(...).
 */
void FluxSurfaceAveragerStellarator::Rebuild(){
    this->nr = rGrid->GetNr();

    // if using fixed quadrature, store all quantities on the theta grid
    if(!integrateAdaptive){
        BOverBmin->InterpolateMagneticDataToThetaPhi(theta, ntheta_interp, phi, nphi_interp);
        Jacobian->InterpolateMagneticDataToThetaPhi(theta, ntheta_interp, phi, nphi_interp);
        BdotGradphi->InterpolateMagneticDataToThetaPhi(theta, ntheta_interp, phi, nphi_interp);
        gttOverJ2->InterpolateMagneticDataToThetaPhi(theta, ntheta_interp, phi, nphi_interp);
        gtpOverJ2->InterpolateMagneticDataToThetaPhi(theta, ntheta_interp, phi, nphi_interp);
        ROverR0->InterpolateMagneticDataToTheta(theta, ntheta_interp);
        NablaR2->InterpolateMagneticDataToTheta(theta, ntheta_interp);
    }

    // Calculate flux-surface averaged Jacobian and hand over to RadialGridStellarator.
    real_t *VpVol   = new real_t[nr];
    real_t *VpVol_f = new real_t[nr+1];    
    for(len_t ir=0; ir<nr;  ir++) 
        VpVol[ir]   = EvaluateFluxSurfaceIntegral(ir, FLUXGRIDTYPE_DISTRIBUTION, RadialGridStellarator::FSA_FUNC_UNITY, nullptr, RadialGridStellarator::FSA_PARAM_UNITY);
    for(len_t ir=0; ir<=nr; ir++)
        VpVol_f[ir] = EvaluateFluxSurfaceIntegral(ir, FLUXGRIDTYPE_RADIAL, RadialGridStellarator::FSA_FUNC_UNITY, nullptr, RadialGridStellarator::FSA_PARAM_UNITY);
    
    rGrid->SetVpVol(VpVol,VpVol_f);
}

/**
 *  Evaluates the flux surface average <F> of a function F = F(B/Bmin, |B \cdot \nabla\varphi|, g_{\theta\theta}/J^2, g_{\theta\theta}/J^2) on radial grid point ir. 
 */
real_t FluxSurfaceAveragerStellarator::CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par, const int_t *F_list){
    real_t VpVol = GetVpVol(ir,fluxGridType);

    // treat singular point r=0 separately where orbit parameters are constant 
    if(VpVol == 0) 
        return F(1,1,1,1,par); // TODO: should we keep this as 1?

    // otherwise use regular method
    return EvaluateFluxSurfaceIntegral(ir,fluxGridType, F, par, F_list) / VpVol;
}

/**
 * The function returns the integrand of the flux surface integral at theta and phi,
 * and is used with the adaptive quadrature
 */
struct FluxSurfaceIntegralParamsThetaPhi {
    real_t(*Function)(real_t,real_t,real_t,real_t,void*); void *FSApar; const int_t *F_list; len_t ir; real_t Bmin; real_t phi;
    FluxSurfaceAveragerStellarator *FSA; fluxGridType fgType;
};
real_t FluxSurfaceAveragerStellarator::FluxSurfaceIntegralFunctionThetaPhi(real_t theta, void *p){
    struct FluxSurfaceIntegralParamsThetaPhi *params = (struct FluxSurfaceIntegralParamsThetaPhi *) p;
    len_t ir = params->ir;
    fluxGridType fluxGridType = params->fgType;
    real_t Bmin = params->Bmin;
    real_t phi = params->phi;

    real_t B,Jacobian,BdotGradphi,gttOverJ2,gtpOverJ2;
    params->FSA->GeometricQuantitiesAtThetaPhi(ir,theta,phi,B,Jacobian,BdotGradphi,gttOverJ2,gtpOverJ2,fluxGridType);
    real_t BOverBmin=1;
    if(Bmin != 0)
        BOverBmin = B/Bmin;
    const int_t *Flist = params->F_list;
    real_t Function = (Flist!=nullptr) ?
        AssembleFSAFunc(BOverBmin, BdotGradphi, gttOverJ2, gtpOverJ2, Flist) 
        : params->Function(BOverBmin, BdotGradphi, gttOverJ2, gtpOverJ2, params->FSApar);  
    return Jacobian*Function;
}

/**
 * The function returns the integrand of the flux surface integral at phi integrated over theta,
 * and is used with the adaptive quadrature
 */
struct FluxSurfaceIntegralParamsPhi {
    real_t(*Function)(real_t,real_t,real_t,real_t,void*); void *FSApar; const int_t *F_list; len_t ir; real_t Bmin;
    FluxSurfaceAveragerStellarator *FSA; fluxGridType fgType;
};
real_t FluxSurfaceAveragerStellarator::FluxSurfaceIntegralFunctionPhi(real_t phi, void *p){
    real_t fluxSurfaceIntegral;
    gsl_function GSL_func; 
    struct FluxSurfaceIntegralParamsPhi *params_phi = (struct FluxSurfaceIntegralParamsPhi *) p;
    FluxSurfaceIntegralParamsThetaPhi params = {params_phi->Function, params_phi->FSApar, params_phi->F_list, params_phi->ir, params_phi->Bmin, phi, params_phi->FSA, params_phi->fgType}; 

    GSL_func.function = &(FluxSurfaceIntegralFunctionThetaPhi);
    GSL_func.params = &params;
    real_t epsabs = 0, epsrel = 1e-4, lim = params_phi->FSA->gsl_adaptive_theta->limit, error;
    gsl_integration_qag(&GSL_func, 0, params_phi->FSA->theta_max,epsabs,epsrel,lim,params_phi->FSA->QAG_KEY,params_phi->FSA->gsl_adaptive_theta,&fluxSurfaceIntegral, &error);
    
    return fluxSurfaceIntegral;
}

real_t FluxSurfaceAveragerStellarator::EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par, const int_t *F_list){
    real_t fluxSurfaceIntegral = 0;

    // Calculate using fixed quadrature:
    if(!integrateAdaptive){
        const real_t *BOverBmin = this->BOverBmin->GetData(ir, fluxGridType);
        const real_t *Jacobian  = this->Jacobian->GetData(ir,fluxGridType);
        const real_t *BdotGradphi = this->BdotGradphi->GetData(ir, fluxGridType);
        const real_t *gttOverJ2   = this->gttOverJ2->GetData(ir, fluxGridType);
        const real_t *gtpOverJ2   = this->gtpOverJ2->GetData(ir, fluxGridType);

        bool hasFlist = (F_list!=nullptr);    
        for (len_t it = 0; it<ntheta_interp; it++)
            for (len_t ip = 0; ip<nphi_interp; ip++)
                fluxSurfaceIntegral += weights_phi[ip] * weights_theta[it] * Jacobian[it*nphi_interp+ip] 
                    * (hasFlist ? AssembleFSAFunc(BOverBmin[it*nphi_interp+ip], BdotGradphi[it*nphi_interp+ip], gttOverJ2[it*nphi_interp+ip], gtpOverJ2[it*nphi_interp+ip], F_list) 
                        : F(BOverBmin[it*nphi_interp+ip], BdotGradphi[it*nphi_interp+ip], gttOverJ2[it*nphi_interp+ip], gtpOverJ2[it*nphi_interp+ip], par));
    // or by using adaptive quadrature:
    } else {
        gsl_function GSL_func; 
        FluxSurfaceIntegralParamsPhi params = {F, par, F_list, ir, GetBmin(ir,fluxGridType), this, fluxGridType}; 
        GSL_func.function = &(FluxSurfaceIntegralFunctionPhi);
        GSL_func.params = &params;
        real_t epsabs = 0, epsrel = 1e-4, lim = gsl_adaptive_phi->limit, error;
        gsl_integration_qag(&GSL_func, 0, phi_max,epsabs,epsrel,lim,QAG_KEY,gsl_adaptive_phi,&fluxSurfaceIntegral, &error);
    }
    return fluxSurfaceIntegral;  
} 



/**
 * Deallocate quadrature.
 */
void FluxSurfaceAveragerStellarator::DeallocateQuadrature(){
    // TODO: Take back if BA
    //gsl_integration_workspace_free(gsl_adaptive);
    //gsl_integration_workspace_free(gsl_adaptive_outer);
    //gsl_integration_qaws_table_free(qaws_table);
    gsl_integration_workspace_free(gsl_adaptive_theta);
    gsl_integration_workspace_free(gsl_adaptive_phi);
    if(gsl_w_theta != nullptr)
        gsl_integration_fixed_free(gsl_w_theta);
    if(gsl_w_phi != nullptr)
        gsl_integration_fixed_free(gsl_w_phi);
}

/*
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
void FluxSurfaceAveragerStellarator::InitializeQuadrature(quadrature_method q_method){
    // TODO: Take back if BA
    //gsl_adaptive = gsl_integration_workspace_alloc(1000);
    //gsl_adaptive_outer = gsl_integration_workspace_alloc(1000);
    gsl_adaptive_theta = gsl_integration_workspace_alloc(1000);
    gsl_adaptive_phi   = gsl_integration_workspace_alloc(1000);
    std::function<real_t(real_t,real_t,real_t)>  QuadWeightFunction;
    theta_max = 2 * M_PI; 
    if(nfp > 0)
        phi_max = M_PI / nfp;
    else 
        phi_max = 2 * M_PI;
    
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
            throw FVMException("Quadrature rule '%d' not supported by FluxSurfaceAveragerStellarator.", q_method);            
    }

    
    gsl_w_theta = gsl_integration_fixed_alloc(quadratureRule,ntheta_interp,0,theta_max,0,0);
    this->theta = gsl_w_theta->x;
    this->weights_theta = gsl_w_theta->weights;

    gsl_w_phi = gsl_integration_fixed_alloc(quadratureRule,nphi_interp,0,phi_max,0,0);
    this->phi = gsl_w_phi->x;
    this->weights_phi = gsl_w_phi->weights;

    // We divide the weights by the QuadWeightFunction to cast the
    // equation on the form required by the quadrature rule.
    for(len_t it=0; it<ntheta_interp; it++)
        weights_theta[it]/= QuadWeightFunction(theta[it],0,theta_max);
    for(len_t ip=0; ip<nphi_interp; ip++)
        weights_phi[ip]/= QuadWeightFunction(phi[ip],0,phi_max);

    // If symmetric field we integrate from 0 to pi and multiply result by 2. 
    if(nfp > 0)
        for(len_t it=0; it<ntheta_interp; it++)
            weights_phi[it] *= 2 * nfp;
}

/** TODO: Take back code for BA, see original FluxSurfaceAverager */

/**
 * Returns the function to be bounce averaged, evaluated at poloidal angle theta.
 * If Flist is provided, returns the function
 *      BA_Func = Flist[5] * xiOverXi0^Flist[0] * BOverBmin^Flist[1] 
 *               * BdotGradphi^Flist[2] * NablaR2^Flist[3],
 */
real_t FluxSurfaceAverager::AssembleBAFunc(real_t xiOverXi0,real_t BOverBmin, real_t BdotGradphi, real_t gttOverJ2, real_t gtpOverJ2, const int_t *Flist){
    real_t BA_Func = Flist[5];
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
            BA_Func *= BdotGradphi;
    else if(Flist[2]<0)
        for(int_t k=0; k<-Flist[2]; k++)
            BA_Func /= BdotGradphi;
    if(Flist[3]>0)
        for(int_t k=0; k<Flist[3]; k++)
            BA_Func *= gttOverJ2;
    else if(Flist[3]<0)
        for(int_t k=0; k<-Flist[3]; k++)
            BA_Func /= gttOverJ2;
    if(Flist[4]>0)
        for(int_t k=0; k<Flist[4]; k++)
            BA_Func *= gtpOverJ2;
    else if(Flist[4]<0)
        for(int_t k=0; k<-Flist[4]; k++)
            BA_Func /= gtpOverJ2;
    return BA_Func;
}

/**
 * Returns the function to be flux surface averaged, evaluated at poloidal angle theta.
 * If Flist is provided, returns the function
 *      FSA_Func = Flist[4] * BOverBmin^Flist[0] 
 *               * BdotGradphi^Flist[1] * gttOverJ2^Flist[2] * gttOverJ2^Flist[3],
 */
real_t FluxSurfaceAveragerStellarator::AssembleFSAFunc(real_t BOverBmin, real_t BdotGradphi, real_t gttOverJ2, real_t gtpOverJ2, const int_t *Flist){
    real_t FSA_Func = Flist[3];
    if(Flist[0]>0)
        for(int_t k=0; k<Flist[0]; k++)
            FSA_Func *= BOverBmin;
    else if(Flist[0]<0)
        for(int_t k=0; k<-Flist[0]; k++)
            FSA_Func /= BOverBmin;
    if(Flist[1]>0)
        for(int_t k=0; k<Flist[1]; k++)
            FSA_Func *= BdotGradphi;
    else if(Flist[1]<0)
        for(int_t k=0; k<-Flist[1]; k++)
            FSA_Func /= BdotGradphi;
    if(Flist[2]>0)
        for(int_t k=0; k<Flist[2]; k++)
            FSA_Func *= gttOverJ2;
    else if(Flist[2]<0)
        for(int_t k=0; k<-Flist[2]; k++)
            FSA_Func /= gttOverJ2;
    if(Flist[3]>0)
        for(int_t k=0; k<Flist[3]; k++)
            FSA_Func *= gtpOverJ2;
    else if(Flist[3]<0)
        for(int_t k=0; k<-Flist[3]; k++)
            FSA_Func /= gtpOverJ2;
    return FSA_Func;
}

/** TODO: This is a bit wrong.... We write the first element twice, but it is never used.
 * Print the variation of B(theta) to stdout.
 */
/*void FluxSurfaceAveragerStellarator::PrintBOfThetaPhi(const len_t ir, const len_t N, enum fluxGridType fgt) {
	printf("B(theta,phi) at ir = " LEN_T_PRINTF_FMT "\n", ir);
	printf("theta = [%.12e", -M_PI);
	for (len_t i = 1; i < N; i++)
		printf(",%.12e", 2*M_PI * (i/((real_t)N)) - M_PI);
	printf("]\n");
	printf("phi = [%.12e", 0.);
	for (len_t i = 1; i < N; i++)
		printf(",%.12e", 2*M_PI * (i/((real_t)N)));
	printf("]\n");
		
	printf("B     = [%.12e", BAtThetaPhi(ir, -M_PI, 0,  fgt));
	for (len_t i = 0; i < N; i++)
        for (len_t j = 0; j < N; j++)
		    printf(",%.12e", BAtThetaPhi(ir, 2*M_PI * (j/((real_t)N)) - M_PI, 2*M_PI * (i/((real_t)N)), fgt));
	printf("]\n");
}*/

