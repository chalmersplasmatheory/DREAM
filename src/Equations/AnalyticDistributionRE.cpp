/**
 * Functions to evaluate the analytic pitch-angle distribution. 
 * (maybe other analytic distributions will be added later?)
 */

#include "DREAM/Equations/AnalyticDistributionRE.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
AnalyticDistributionRE::AnalyticDistributionRE(
    FVM::RadialGrid *rGrid, FVM::UnknownQuantityHandler *u, PitchScatterFrequency *nuD, 
    CollisionQuantity::collqty_settings *cqset, dist_mode mode, 
    real_t thresholdToNeglectTrappedContribution
) : AnalyticDistribution(rGrid, u), nuD(nuD), collSettings(cqset), mode(mode), thresholdToNeglectTrappedContribution(thresholdToNeglectTrappedContribution){
    this->gsl_ad_w = gsl_integration_workspace_alloc(1000);
    this->id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
}

/**
 * Destructor
 */
AnalyticDistributionRE::~AnalyticDistributionRE(){
    gsl_integration_workspace_free(gsl_ad_w);
}

/**
 * Same as evaluatePitchDistribution but takes A (width parameter) instead of p, E 
 * and inSettings used to create look-up-table in the Eceff calculation.
 */
real_t AnalyticDistributionRE::evaluatePitchDistributionFromA(
    len_t ir, real_t xi0, real_t A
){
    if(mode == RE_PITCH_DIST_SIMPLE)
        return evaluateApproximatePitchDistributionFromA(ir,xi0,A);
    else
        return evaluateAnalyticPitchDistributionFromA(ir,xi0,A);
}

/**
 * Returns xi0/<xi> (the integral of which appears in AnalyticPitchDistribution).
 */
struct distExponentParams {len_t ir; FVM::RadialGrid *rGrid;};
real_t distExponentIntegral(real_t xi0, void *par){
    struct distExponentParams *params = (struct distExponentParams *) par;
    std::function<real_t(real_t,real_t,real_t)> xiFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return sqrt(1 - BOverBmin*(1-xi0*xi0));};
    len_t ir = params->ir;
    FVM::RadialGrid *rGrid = params->rGrid;
    real_t signXi0 = ( (xi0>0) - (xi0<0));
    real_t xiAvg = signXi0 * rGrid->CalculateFluxSurfaceAverage(ir,FVM::FLUXGRIDTYPE_DISTRIBUTION, xiFunc);
    return xi0/xiAvg;
}

/**
 * Calculates the (semi-)analytic pitch-angle distribution predicted in the 
 * near-threshold regime, where the momentum flux is small compared 
 * to the characteristic pitch flux, and we obtain the approximate 
 * kinetic equation phi_xi = 0.
 */
real_t AnalyticDistributionRE::evaluateAnalyticPitchDistributionFromA(
    len_t ir, real_t xi0, real_t A
){
    real_t xiT = rGrid->GetXi0TrappedBoundary(ir);  

    // This block carries defines the integration int(xi0/<xi(xi0)> dxi0, xi1, x2) 
    //////////////////////////////
    gsl_function GSL_func;
    distExponentParams params = {ir,rGrid};
    GSL_func.function = &(distExponentIntegral);
    GSL_func.params = &params;
    real_t abserr;
    real_t epsabs = 0, epsrel = 3e-3, lim = gsl_ad_w->limit;
    #define F(xi1,xi2,pitchDist) gsl_integration_qags(&GSL_func, xi1,xi2,epsabs,epsrel,lim,gsl_ad_w, &pitchDist, &abserr)
    //////////////////////////////    

    real_t dist1 = 0;
    real_t dist2 = 0;

    if ( (xi0>xiT) || (xiT<thresholdToNeglectTrappedContribution) )
        F(xi0,1.0,dist1);
    else if ( (-xiT <= xi0) && (xi0 <= xiT) )
        F(xiT,1.0,dist1);
    else{ // (xi0 < -xiT)
        F(xi0,-xiT,dist1);
        F(xiT,1.0,dist2);
    }
    
    #undef F
    
    return exp(-A*(dist1+dist2));
}

/**
 * Same as evaluteAnalyticPitchDistribution, but approximating
 * xi0/<xi> = 1 for passing and 0 for trapped (thus avoiding the 
 * need for the numerical integration).
 */
real_t AnalyticDistributionRE::evaluateApproximatePitchDistributionFromA(len_t ir, real_t xi0, real_t A){
    real_t xiT = rGrid->GetXi0TrappedBoundary(ir);
    real_t dist1 = 0;
    real_t dist2 = 0;

    if ( (xi0>xiT) || (xiT<this->thresholdToNeglectTrappedContribution) )
        dist1 = 1-xi0;
    else if ( (-xiT <= xi0) && (xi0 <= xiT) )
        dist1 = 1-xiT;
    else{ // (xi0 < -xiT)
        dist1 = 1-xiT;
        dist2 = -xiT - xi0;
    }
    return exp(-A*(dist1+dist2));
}

//                                                     (len_t ir, real_t xi0, real_t p, real_t *dfdxi0, real_t *dfdp, real_t *dfdr)
real_t AnalyticDistributionRE::evaluateFullDistribution(len_t   , real_t    , real_t  , real_t *      , real_t *    , real_t *){
    return NAN;
} 

//                                                       (len_t ir, real_t p, real_t *dfdp, real_t *dfdr)
real_t AnalyticDistributionRE::evaluateEnergyDistribution(len_t,    real_t ,  real_t *,     real_t *){
    return NAN;
    // implement avalanche distribution
}

real_t AnalyticDistributionRE::evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t *dfdxi0, real_t *dfdp, real_t *dfdr){
    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    real_t Eterm = unknowns->GetUnknownData(id_Eterm)[ir];
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 
    real_t pNuD = p*nuD->evaluateAtP(ir,p,collSettings);    
    real_t A = 2*E/pNuD;

    if(dfdxi0!=nullptr){
        // evaluate pitch derivative
    }
    if(dfdp!=nullptr){
        //evaluate p derivative
    }
    if(dfdr!=nullptr){
        //evaluate r derivative
    }
    return evaluatePitchDistributionFromA(ir, xi0, A);
}
