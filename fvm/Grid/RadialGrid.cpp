/**
 * Implementation of the radial grid.
 */

#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "gsl/gsl_integration.h"

using namespace std;
using namespace DREAM::FVM;

/***********************
 * Constructors        *
 ***********************/
/**
 * Initialize an empty grid by only specifying the
 * grid size.
 *
 * rg: Object to use for (re-)generating the radial grid.
 * t0: Time to initialize grid at.
 * ntheta_interp: Poloidal angle resolution in quadrature for flux surface and bounce averages.
 */
RadialGrid::RadialGrid(RadialGridGenerator *rg, const real_t /*t0*/,
    FluxSurfaceAverager::interp_method im, FluxSurfaceAverager::quadrature_method qm_passing)
    : nr(rg->GetNr()), generator(rg) {

    bool geometryIsSymmetric = rg->IsFieldSymmetric();
    len_t ntheta_interp_passing = rg->GetNthetaInterp();

    fluxSurfaceAverager = new FluxSurfaceAverager(this,rg,geometryIsSymmetric,ntheta_interp_passing,im,qm_passing);
}

/**
 * Destructor.
 */
RadialGrid::~RadialGrid(){    
    DeallocateGrid();
    DeallocateFSAvg();
    
    DeallocateReferenceMagneticData();
    DeallocateMagneticExtremumData();
    
    if(this->VpVol!=nullptr){
        delete [] this->VpVol;
        delete [] this->VpVol_f;
    }
	if (this->psiToroidal != nullptr) {
		delete [] this->psiToroidal;
		delete [] this->psiToroidal_f;
	}

    delete this->generator;
    delete this->fluxSurfaceAverager;
}

/***************************
 * PUBLIC METHODS          *
 ***************************/
/**
 * Rebuilds any non-static (i.e. time dependent) grids
 * used. This can be used if, for example, a dynamically
 * evolving magnetic equilibrium is used, or if some
 * grids are adaptive.
 *
 * t: Time to which re-build the grids for.
 */
bool RadialGrid::Rebuild(const real_t t) {
    // Re-build radial grid
    if (this->generator->NeedsRebuild(t)){
        bool rebuilt = this->generator->Rebuild(t, this);
        return rebuilt;
    }
    else return false;
}

/**
 * Rebuilds magnetic-field data and initializes flux 
 * surface average calculations.
 */
void RadialGrid::RebuildJacobians(){ 
    this->generator->RebuildJacobians(this);
    fluxSurfaceAverager->Rebuild();
    RebuildFluxSurfaceAveragedQuantities();
}


/**
 * Grid deallocator
 */
void RadialGrid::DeallocateGrid() {
    if (this->r == nullptr)
        return;

    delete [] this->dr_f;
    delete [] this->dr;
    delete [] this->r_f;
    delete [] this->r;
}


/**
 * Calculate flux surface average
 */
real_t RadialGrid::CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,void*), void *par, const int_t *Flist){
    return fluxSurfaceAverager->CalculateFluxSurfaceAverage(ir, fluxGridType, F, par, Flist);
}
/**
 * Evaluate flux surface integral
 */
real_t RadialGrid::EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,void*), void *par, const int_t *Flist){
    return fluxSurfaceAverager->EvaluateFluxSurfaceIntegral(ir, fluxGridType, F, par, Flist);
}
/**
 * Calculate bounce average at arbitrary p and xi
 */
real_t RadialGrid::CalculatePXiBounceAverageAtP(len_t ir, real_t xi0, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par, const int_t *Flist){
    return fluxSurfaceAverager->CalculatePXiBounceAverageAtP(ir,xi0,fluxGridType,F, par, Flist);
}
/**
 * Evaluate bounce integral at arbitrary p and xi
 */
real_t RadialGrid::EvaluatePXiBounceIntegralAtP(len_t ir, real_t xi0, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par, const int_t *Flist){
    return fluxSurfaceAverager->EvaluatePXiBounceIntegralAtP(ir,xi0,fluxGridType,F, par, Flist);
}



/**
 * Sets magnetic field quantities that have been
 * generated by a RadialGridGenerator (FluxSurfaceAverager
 * takes ownership of some and will deallocate them).
 */
void RadialGrid::SetMagneticExtremumData(
    real_t *Bmin, real_t *Bmin_f,
    real_t *Bmax, real_t *Bmax_f,
    real_t *theta_Bmin, real_t *theta_Bmin_f,
    real_t *theta_Bmax, real_t *theta_Bmax_f,
    real_t *xi0TrappedBoundary, real_t *xi0TrappedBoundary_f
){
    DeallocateMagneticExtremumData();

    this->Bmin           = Bmin;
    this->Bmin_f         = Bmin_f;
    this->Bmax           = Bmax;
    this->Bmax_f         = Bmax_f;
    this->xi0TrappedBoundary = xi0TrappedBoundary;
    this->xi0TrappedBoundary_f = xi0TrappedBoundary_f;
    
    fluxSurfaceAverager->SetReferenceMagneticFieldData(
        theta_Bmin, theta_Bmin_f,
        theta_Bmax, theta_Bmax_f
    );
}

/**
 * Sets reference magnetic field quantities that have been
 * generated by a RadialGridGenerator.
 */
void RadialGrid::SetReferenceMagneticFieldData(
    real_t *BtorGOverR0, real_t *BtorGOverR0_f,
    real_t *psiPrimeRef, real_t *psiPrimeRef_f, 
    real_t R0
){
    DeallocateReferenceMagneticData();

    this->BtorGOverR0    = BtorGOverR0;
    this->BtorGOverR0_f  = BtorGOverR0_f;
    this->psiPrimeRef    = psiPrimeRef;
    this->psiPrimeRef_f  = psiPrimeRef_f;
    this->R0             = R0;
}


/**
 * Calculate and store flux surface averages.
 */
void RadialGrid::RebuildFluxSurfaceAveragedQuantities(){
 real_t 
    *effectivePassingFraction   = nullptr, 
    *effectivePassingFraction_f = nullptr, 
    *FSA_B2   = nullptr,
    *FSA_B2_f = nullptr,
    *FSA_B    = nullptr,
    *FSA_B_f  = nullptr,
    *FSA_nablaR2OverR2   = nullptr,
    *FSA_nablaR2OverR2_f = nullptr, 
    *FSA_1OverR2   = nullptr,
    *FSA_1OverR2_f = nullptr;

    SetFluxSurfaceAverage(FSA_1OverR2,FSA_1OverR2_f, FSA_FUNC_ONE_OVER_R_SQUARED, nullptr, FSA_PARAM_ONE_OVER_R_SQUARED);
    SetFluxSurfaceAverage(FSA_B,FSA_B_f, FSA_FUNC_B, nullptr, FSA_PARAM_B);
    SetFluxSurfaceAverage(FSA_B2,FSA_B2_f, FSA_FUNC_B_SQUARED, nullptr, FSA_PARAM_B_SQUARED);
    SetFluxSurfaceAverage(FSA_nablaR2OverR2,FSA_nablaR2OverR2_f, FSA_FUNC_NABLA_R_SQUARED_OVER_R_SQUARED, nullptr, FSA_PARAM_NABLA_R_SQUARED_OVER_R_SQUARED);
    
    SetEffectivePassingFraction(effectivePassingFraction,effectivePassingFraction_f, FSA_B2, FSA_B2_f);

    InitializeFSAvg(effectivePassingFraction,effectivePassingFraction_f,
        FSA_B,FSA_B_f,FSA_B2,FSA_B2_f,FSA_1OverR2, FSA_1OverR2_f,FSA_nablaR2OverR2,FSA_nablaR2OverR2_f);

    // set toroidal flux psi_t defined by dpsi_t/dpsi_p = qR0 (safety factor)
    // or equivalently as the toroidal magnetic field integrated over a 
    // poloidal cross section
    if(this->psiToroidal != nullptr) {
        delete [] psiToroidal;
        delete [] psiToroidal_f;
    }

    this->psiToroidal = new real_t[nr];
    this->psiToroidal_f = new real_t[nr+1];

    // Calculate psi_t by integrating using a trapezoidal rule.
    real_t x, x_f = (1. /(2. * M_PI)) * VpVol_f[0]*BtorGOverR0_f[0]*FSA_1OverR2_f[0];
    psiToroidal_f[0] = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        x    = (1. /(2. * M_PI)) * VpVol[ir]*BtorGOverR0[ir]*FSA_1OverR2[ir];
        psiToroidal[ir] = psiToroidal_f[ir] + 0.5*(x_f+x)*(r[ir]-r_f[ir]);

        x_f  = (1. /(2. * M_PI)) * VpVol_f[ir+1]*BtorGOverR0_f[ir+1]*FSA_1OverR2_f[ir+1];
        psiToroidal_f[ir+1] = psiToroidal[ir] + 0.5*(x_f+x)*(r_f[ir+1]-r[ir]);
    }
}

/**
 * Helper method to store flux surface averages.
 */
void RadialGrid::SetFluxSurfaceAverage(real_t *&FSA_quantity, real_t *&FSA_quantity_f, real_t(*F)(real_t,real_t,real_t,void*), void *par, const int_t *Flist){
    FSA_quantity   = new real_t[GetNr()];
    FSA_quantity_f = new real_t[GetNr()+1];

    for(len_t ir=0; ir<nr; ir++)
        FSA_quantity[ir] = CalculateFluxSurfaceAverage(ir, FLUXGRIDTYPE_DISTRIBUTION, F, par, Flist);

    for(len_t ir=0; ir<=nr; ir++)
        FSA_quantity_f[ir] = CalculateFluxSurfaceAverage(ir, FLUXGRIDTYPE_RADIAL, F, par, Flist);
}


/**
 * The function is used in the evaluation of the effective passing fraction,
 * and represents x / <1-x B/Bmax>
 */
real_t RadialGrid::effectivePassingFractionIntegrand(real_t x, void *p){
    struct EPF_params *params = (struct EPF_params *) p;
    params->x = x;
    RadialGrid *rGrid = params->rGrid;
    len_t ir = params->ir;
    fluxGridType fluxGridType = params->fgType;
    return x/ rGrid->CalculateFluxSurfaceAverage(ir, fluxGridType, FSA_FUNC_EFF_PASS_FRAC, params);
}

/**
 * Calculates and stores the effective fraction of passing electrons
 */
void RadialGrid::SetEffectivePassingFraction(real_t *&EPF, real_t *&, real_t *FSA_B2, real_t*){
    gsl_integration_workspace *gsl_w = gsl_integration_workspace_alloc(1000);
    gsl_function EPF_func;
    EPF_params paramstruct;
    real_t EPF_integral;
    EPF = new real_t[GetNr()];
    real_t error = 0;
    for (len_t ir=0; ir<GetNr(); ir++){
        real_t Bmin = GetBmin(ir);
        real_t Bmax = GetBmax(ir);
        real_t BminOverBmax;
        if(Bmin==Bmax) // handles B(theta) = 0 case
            BminOverBmax = 1; 
        else
            BminOverBmax = Bmin/Bmax;
        paramstruct = {-1,BminOverBmax,ir,this,FLUXGRIDTYPE_DISTRIBUTION};
        EPF_func.function = &(effectivePassingFractionIntegrand);
        EPF_func.params = &paramstruct;
        real_t epsabs = 0, epsrel = 1e-4, lim = gsl_w->limit;
        gsl_integration_qags(&EPF_func, 0,1,epsabs,epsrel,lim,gsl_w,&EPF_integral, &error);
        EPF[ir] = (3.0/4.0) * BminOverBmax*BminOverBmax*FSA_B2[ir] * EPF_integral;
    }
    gsl_integration_workspace_free(gsl_w);
}

/**
 * Set flux surface averages
 */
void RadialGrid::InitializeFSAvg(
    real_t *epf, real_t *epf_f, real_t *Bavg, real_t *Bavg_f, 
    real_t *B2avg, real_t *B2avg_f,
    real_t *OneOverR2_avg, real_t *OneOverR2_avg_f,
    real_t *nablaR2OverR2_avg, real_t *nablaR2OverR2_avg_f
){
    DeallocateFSAvg();
    this->effectivePassingFraction   = epf;
    this->effectivePassingFraction_f = epf_f;
    this->FSA_B                      = Bavg;
    this->FSA_B_f                    = Bavg_f;
    this->FSA_B2                     = B2avg;
    this->FSA_B2_f                   = B2avg_f;
    this->FSA_1OverR2                = OneOverR2_avg;
    this->FSA_1OverR2_f              = OneOverR2_avg_f;
    this->FSA_nablaR2OverR2          = nablaR2OverR2_avg;
    this->FSA_nablaR2OverR2_f        = nablaR2OverR2_avg_f;   
}


real_t RadialGrid::GetFluxSurfaceR(len_t ir) {
    return this->GetR0() + this->GetGridGenerator()->GetFluxSurfaceRMinusR0()[ir];
}

real_t RadialGrid::GetFluxSurfaceZ(len_t ir) {
    return this->GetGridGenerator()->GetFluxSurfaceZMinusZ0()[ir];
}



/**
 * Deallocate flux surface averages
 */
void RadialGrid::DeallocateFSAvg(){
    if (this->effectivePassingFraction == nullptr)
        return;

    delete [] this->FSA_B;
    delete [] this->FSA_B_f;
    delete [] this->FSA_B2;
    delete [] this->FSA_B2_f;
    delete [] this->effectivePassingFraction;
    delete [] this->FSA_nablaR2OverR2;
    delete [] this->FSA_nablaR2OverR2_f;
    delete [] this->FSA_1OverR2;
    delete [] this->FSA_1OverR2_f;
}

