/**
 * Functions to evaluate the analytic pitch-angle distribution. 
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
) : AnalyticDistribution(rGrid, u), nuD(nuD), collSettings(cqset), mode(mode), 
    thresholdToNeglectTrappedContribution(thresholdToNeglectTrappedContribution),
    id_Eterm(u->GetUnknownID(OptionConstants::UQTY_E_FIELD))
{
    this->gsl_ad_w = gsl_integration_workspace_alloc(1000);

    GridRebuilt();
}


void AnalyticDistributionRE::Deallocate(){
    if(xi0OverXiSpline != nullptr){
        for(len_t ir=0; ir<nr; ir++){
            gsl_spline_free(xi0OverXiSpline[ir]);
            gsl_interp_accel_free(xiSplineAcc[ir]);
        }
        delete [] xi0OverXiSpline;
        delete [] xiSplineAcc;
        delete [] integralOverFullPassing;
    }
    if(REDistNormFactor_Spline != nullptr){
        for(len_t ir=0; ir<nr; ir++){
            gsl_spline_free(REDistNormFactor_Spline[ir]);
            gsl_interp_accel_free(REDistNormFactor_Accel[ir]);
        }
        delete [] REDistNormFactor_Spline;
        delete [] REDistNormFactor_Accel;
    }

    
}
/**
 * Destructor
 */
AnalyticDistributionRE::~AnalyticDistributionRE(){
    Deallocate();
    gsl_integration_workspace_free(gsl_ad_w);
}

bool AnalyticDistributionRE::GridRebuilt(){
    Deallocate();
    this->AnalyticDistribution::GridRebuilt();

    if(mode==RE_PITCH_DIST_FULL)
        constructXiSpline();

    constructVpSplines();

    return true;
}

/**
 * Sets splines for the pitch-dist-bounce-integrated metric
 *      VpRE = int( Vp * f(xi) dxi0)
 * where f is the 'Lehtinen theory' pitch distribution.
 * VpRE is therefore also the normalization factor such that
 * the distribution f/VpRE integrates over xi0 to unity, and
 * the jacobian for the remaining momentum distribution is 2*pi*p^2VpVol. 
 */
void AnalyticDistributionRE::constructVpSplines(){
    std::function<real_t(real_t)> trivialFunc = [](real_t){return 1.0;};
    const len_t N_VP_SPLINE = 100;
    const len_t N_RE_DIST_SPLINE = 50;

    real_t 
        *xArray = new real_t[N_RE_DIST_SPLINE],
        *REDistIntegralArray = new real_t[N_RE_DIST_SPLINE];

    REPitchDistributionAveragedBACoeff::GenerateNonUniformXArray(xArray, N_RE_DIST_SPLINE);

    REDistNormFactor_Accel  = new gsl_interp_accel*[nr];
    REDistNormFactor_Spline = new gsl_spline*[nr];
    for(len_t ir=0; ir<nr; ir++){
        real_t xiT = rGrid->GetXi0TrappedBoundary(ir);
        // Create a temporary spline over Vp:
        gsl_spline *VpSpline;
        gsl_interp_accel *VpAcc = gsl_interp_accel_alloc();;
        REPitchDistributionAveragedBACoeff::GenerateBASpline(ir, rGrid,
            xiT, N_VP_SPLINE, FVM::RadialGrid::BA_FUNC_UNITY, nullptr, 
            FVM::RadialGrid::BA_PARAM_UNITY, VpSpline
        );

        // test that VpSpline correctly integrates to 4*pi*VpVol
        /*
        real_t splineIntegral = gsl_spline_eval_integ(VpSpline,0.0,xiT,VpAcc) + 2*gsl_spline_eval_integ(VpSpline,xiT,1.0,VpAcc);
        real_t splineError = fabs(splineIntegral / (4*M_PI*rGrid->GetVpVol(ir)) - 1.0);
        if(splineError > 0.01)
            printf("AnalyticDistributionRE: The integrated VpSpline incurs an error: %f.\n", splineError);
        */

        // Generate spline in A of the normalization factor
        REDistNormFactor_Accel[ir] = gsl_interp_accel_alloc();
        REDistNormFactor_Spline[ir] = gsl_spline_alloc(gsl_interp_steffen, N_RE_DIST_SPLINE);
        REPitchDistributionAveragedBACoeff::ParametersForREPitchDistributionIntegral params 
            = {ir, xiT, 0, this, VpSpline, VpAcc, trivialFunc};
        for(len_t i=0; i<N_RE_DIST_SPLINE; i++){
            real_t A = REPitchDistributionAveragedBACoeff::GetAFromX(xArray[i]);
            params.A = A;
            REDistIntegralArray[i] = REPitchDistributionAveragedBACoeff::EvaluateREDistBounceIntegral(params, gsl_ad_w);
        }
        gsl_spline_init(REDistNormFactor_Spline[ir], xArray, REDistIntegralArray, N_RE_DIST_SPLINE);

        gsl_spline_free(VpSpline);
        gsl_interp_accel_free(VpAcc);
    }
    delete [] xArray;
    delete [] REDistIntegralArray;
}

/**
 * Generates splines over xi0/<xi0> on a xi0 grid 
 */
void AnalyticDistributionRE::constructXiSpline(){
    xi0OverXiSpline = new gsl_spline*[nr];
    xiSplineAcc     = new gsl_interp_accel*[nr];
    real_t *xiArr   = new real_t[N_SPLINE];
    real_t *FuncArr = new real_t[N_SPLINE];
    integralOverFullPassing = new real_t[nr];
    // generate pitch grid for the spline
    for(len_t ir=0; ir<nr; ir++){
        xiSplineAcc[ir]     = gsl_interp_accel_alloc();
        xi0OverXiSpline[ir] = gsl_spline_alloc (gsl_interp_steffen, N_SPLINE);
        real_t xiT  = rGrid->GetXi0TrappedBoundary(ir);
        if(xiT==0) // cylindrical geometry - skip remainder since these splines will not be used
            continue;
        for(len_t k=0; k<N_SPLINE; k++){
            // create uniform xi0 grid on [xiT,1] 
            real_t xi0 = xiT + k*(1.0-xiT)/(N_SPLINE-1);
            xiArr[k]   = xi0;
            // evaluate xi0/<xi> values
            FuncArr[k] = xi0 / rGrid->CalculateFluxSurfaceAverage(
                ir,FVM::FLUXGRIDTYPE_DISTRIBUTION, FVM::RadialGrid::FSA_FUNC_XI, &xi0
            );
        }
        gsl_spline_init (xi0OverXiSpline[ir], xiArr, FuncArr, N_SPLINE);
        // the integral int( xi0/<xi>, xiT, 1 ) over the entire spline 
        // will appear repeatedly and is therefore stored 
        integralOverFullPassing[ir] = gsl_spline_eval_integ(xi0OverXiSpline[ir],xiT,1.0,xiSplineAcc[ir]);
    }
    delete [] xiArr;
    delete [] FuncArr;
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
 * Calculates the (semi-)analytic pitch-angle distribution predicted in the 
 * near-threshold regime, where the momentum flux is small compared 
 * to the characteristic pitch flux, and we obtain the approximate 
 * kinetic equation phi_xi = 0.
 */
real_t AnalyticDistributionRE::evaluateAnalyticPitchDistributionFromA(
    len_t ir, real_t xi0, real_t A
){
    real_t xiT = rGrid->GetXi0TrappedBoundary(ir); 
    if(xiT==0)
        return exp(-A*(1-xi0));

    #define F(xi1,xi2,val) gsl_spline_eval_integ_e(xi0OverXiSpline[ir],xi1,xi2,xiSplineAcc[ir],&val)

    real_t dist1 = 0; // contribution to exponent from positive pitch 
    real_t dist2 = 0; // contribution to exponent from negative pitch

    if (xi0>xiT)
        gsl_spline_eval_integ_e(xi0OverXiSpline[ir],xi0,1.0,xiSplineAcc[ir],&dist1);
    else 
        dist1 = integralOverFullPassing[ir]; // equivalent to F(xiT,1.0,&dist1)
    
    if(xi0<-xiT) // add [xi0, -xiT] part but mirror the interval: spline is 
                 // only defined for positive pitch since it's symmetric
        gsl_spline_eval_integ_e(xi0OverXiSpline[ir],xiT,-xi0,xiSplineAcc[ir],&dist2);
    
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

real_t AnalyticDistributionRE::evaluatePitchDistribution(
    len_t ir, real_t xi0, real_t p, 
    real_t *dfdxi0, real_t *dfdp, real_t *dfdr
){
    real_t A = GetAatP(ir,p);

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

/**
 * Evaluates the pitch distribution width parameter 'A'
 */
real_t AnalyticDistributionRE::GetAatP(len_t ir,real_t p){
    real_t Eterm = unknowns->GetUnknownData(id_Eterm)[ir];
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(rGrid->GetFSA_B2(ir)); 
    real_t pNuD = p*nuD->evaluateAtP(ir,p,collSettings);    
    return 2*E/pNuD;
}

/**
 * Evaluates the "pitch-distribution-bounce jacobian"
 *   VpRE = int((Vp/p^2) * exp(-A*g) dxi0,-1,1)
 */
real_t AnalyticDistributionRE::EvaluateVpREAtA(len_t ir, real_t A){
    return gsl_spline_eval(
        REDistNormFactor_Spline[ir], 
        REPitchDistributionAveragedBACoeff::GetXFromA(A),
        REDistNormFactor_Accel[ir]
    );
}
