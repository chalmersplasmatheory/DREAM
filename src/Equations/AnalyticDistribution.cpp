/**
 * Functions to evaluate the analytic pitch-angle distribution. 
 * (maybe other analytic distributions will be added later?)
 */

#include "DREAM/Equations/AnalyticDistribution.hpp"
#include "FVM/config.h"
#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

using namespace DREAM;

/**
 * Constructor.
 *
 */
AnalyticDistribution::AnalyticDistribution(FVM::RadialGrid *rGrid, PitchScatterFrequency *nuD, OptionConstants::collqty_Eceff_mode Eceff_mode) 
: rGrid(rGrid), nuD(nuD), Eceff_mode(Eceff_mode){}

// ok, maybe re-name the constant so it doen't belong to Eceff?
real_t AnalyticDistribution::evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, 
real_t Eterm, CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w){
    if(Eceff_mode == OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE)
        return evaluateApproximatePitchDistribution(ir,xi0,p,Eterm,inSettings);
    else
        return evaluateAnalyticPitchDistribution(ir,xi0,p,Eterm,inSettings,gsl_ad_w);
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
real_t AnalyticDistribution::evaluateAnalyticPitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, 
CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w){
    const real_t Bmin = rGrid->GetBmin(ir);
    const real_t Bmax = rGrid->GetBmax(ir);
    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 

    real_t pNuD = p*nuD->evaluateAtP(ir,p,inSettings);    
    real_t A = 2*E/pNuD;

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

    real_t thresholdToNeglectTrappedContribution = 1e-6;
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
real_t AnalyticDistribution::evaluateApproximatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings){
    const real_t Bmin = rGrid->GetBmin(ir);
    const real_t Bmax = rGrid->GetBmax(ir);
    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 

    real_t pNuD = p*nuD->evaluateAtP(ir,p,inSettings);    
    real_t A = 2*E/pNuD;

    real_t dist1 = 0;
    real_t dist2 = 0;

    real_t thresholdToNeglectTrappedContribution = 1e-6;
    if ( (xi0>xiT) || (xiT<thresholdToNeglectTrappedContribution) )
        dist1 = 1-xi0;
    else if ( (-xiT <= xi0) && (xi0 <= xiT) )
        dist1 = 1-xiT;
    else{ // (xi0 < -xiT)
        dist1 = 1-xiT;
        dist2 = -xiT - xi0;
    }
        
    return exp(-A*(dist1+dist2));
}
