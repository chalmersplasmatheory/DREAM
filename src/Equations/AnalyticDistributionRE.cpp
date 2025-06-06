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
    id_Eterm(u->GetUnknownID(OptionConstants::UQTY_E_FIELD)),
    id_nre(u->GetUnknownID(OptionConstants::UQTY_N_RE))
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
    std::function<real_t(real_t)> identityFunc = [](real_t){return 1.0;};
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
        REPitchDistributionAveragedBACoeff::GenerateBASpline(
            ir, rGrid, xiT, N_VP_SPLINE, FVM::RadialGrid::BA_FUNC_UNITY, 
            nullptr, FVM::RadialGrid::BA_PARAM_UNITY, VpSpline
        );

        // Generate spline in A of the normalization factor
        REDistNormFactor_Accel[ir] = gsl_interp_accel_alloc();
        REDistNormFactor_Spline[ir] = gsl_spline_alloc(gsl_interp_steffen, N_RE_DIST_SPLINE);
        REPitchDistributionAveragedBACoeff::ParametersForREPitchDistributionIntegral params 
            = {ir, xiT, 0, this, VpSpline, VpAcc, identityFunc};
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
            real_t xi0 = xiT + k*(1.0-xiT)/(N_SPLINE-1.0);
            xiArr[k]   = xi0;
            // evaluate xi0/<xi> values
            FuncArr[k] = xi0 / rGrid->CalculateFluxSurfaceAverage(
                ir,FVM::FLUXGRIDTYPE_DISTRIBUTION, FVM::RadialGrid::FSA_FUNC_XI, &xi0
            );
        }
        xiArr[N_SPLINE-1] = 1.0;
        gsl_spline_init (xi0OverXiSpline[ir], xiArr, FuncArr, N_SPLINE);
        // the integral int( xi0/<xi>, xiT, 1 ) over the entire spline 
        // will appear repeatedly and is therefore stored 
        integralOverFullPassing[ir] = gsl_spline_eval_integ(xi0OverXiSpline[ir],xiArr[0],1.0,xiSplineAcc[ir]);
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
    if(xiT<this->thresholdToNeglectTrappedContribution)
        return exp(-A*(1-xi0));

    real_t dist1 = 0; // contribution to exponent from positive pitch 
    real_t dist2 = 0; // contribution to exponent from negative pitch

    if (xi0>xiT)
        dist1 = gsl_spline_eval_integ(xi0OverXiSpline[ir],xi0,1.0,xiSplineAcc[ir]);
    else 
        dist1 = integralOverFullPassing[ir]; // equivalent to F(xiT,1.0,&dist1)
    
    if(xi0<-xiT) // add [xi0, -xiT] part but mirror the interval: spline is 
                 // only defined for positive pitch since it's symmetric
        dist2 = gsl_spline_eval_integ(xi0OverXiSpline[ir],xiT,-xi0,xiSplineAcc[ir]);
    
    return exp(-A*(dist1+dist2));
}

/**
 * Same as evaluteAnalyticPitchDistribution, but approximating
 * xi0/<xi> = 1 for passing and 0 for trapped (thus avoiding the 
 * need for the numerical integration).
 */
real_t AnalyticDistributionRE::evaluateApproximatePitchDistributionFromA(len_t ir, real_t xi0, real_t A){
    real_t xiT = rGrid->GetXi0TrappedBoundary(ir);
    if(xiT<this->thresholdToNeglectTrappedContribution)
        return exp(-A*(1-xi0));
    real_t dist1 = 0;
    real_t dist2 = 0;

    if(xi0>xiT)
        dist1 = 1-xi0;
    else 
        dist1 = 1-xiT;
    if(xi0<-xiT)
        dist2 = -xiT - xi0;

    return exp(-A*(dist1+dist2));
}

/**
 * Evaluates the energy distribution defined such that 
 *  <n_re> = int(p^2 * [distribution] dp, 0, inf)
 * Assumes that the friction in the acceleration function 'U(p)'
 * is such that it reduces the effective accelerating 
 * electric field by a constant, Eceff, and employs the 
 * zero-pitch-angle limit for the E-field term.
 */
real_t AnalyticDistributionRE::evaluateEnergyDistribution(len_t ir, real_t p, real_t *, real_t *){
    // implement avalanche distribution
    real_t Eterm = unknowns->GetUnknownData(id_Eterm)[ir];
	return evaluateEnergyDistributionWithE(ir, p, Eterm);
}

real_t AnalyticDistributionRE::evaluateEnergyDistributionWithE(
	len_t ir, real_t p, real_t Eterm,
	real_t*, real_t*
) {
    real_t n_re  = unknowns->GetUnknownData(id_nre)[ir];
    real_t Eceff = REFluid->GetEffectiveCriticalField(ir);
    real_t GammaAva = REFluid->GetAvalancheGrowthRate(ir);

    real_t p0 = Constants::ec/(Constants::me*Constants::c) * (fabs(Eterm) - Eceff) * sqrt(rGrid->GetFSA_B2(ir)) / GammaAva;
    real_t F0 = n_re /(p0*p*p) * exp(-p/p0);

    return F0;
}

/**
 * Pitch distribution normalized such that its integral weighted by
 * Vp/(p^2*VpVol) yields unity
 */
real_t AnalyticDistributionRE::evaluatePitchDistribution(
    len_t ir, real_t xi0, real_t p, 
    real_t * /*dfdxi0*/, real_t * /*dfdp*/, real_t * /*dfdr*/
    /*, real_t *dfdA */
){
    real_t A = GetAatP(ir,p);

    real_t D = evaluatePitchDistributionFromA(ir, xi0, A) * rGrid->GetVpVol(ir) / EvaluateVpREAtA(ir, A);
    /* // Implement as need arises
    if(dfdxi0!=nullptr){
        // evaluate pitch derivative
    }
    if(dfdp!=nullptr){
        //evaluate p derivative
    }
    if(dfdr!=nullptr){
        //evaluate r derivative
    }
    if(dfdA!=nullptr){
        //evaluate derivative with respect to A
        real_t h = 1e-3*A;
        real_t Dshift = evaluatePitchDistributionFromA(ir, xi0, A+h) * rGrid->GetVpVol(ir) / EvaluateVpREAtA(ir, A+h);
        *dfdA = (Dshift - D) / h;
    }
    */
    return D;
}

real_t AnalyticDistributionRE::evaluatePitchDistributionWithE(
    len_t ir, real_t xi0, real_t p, real_t E,
    real_t * /*dfdxi0*/, real_t * /*dfdp*/, real_t * /*dfdr*/
    /*, real_t *dfdA */
) {
	real_t A = GetAatP(ir, p, this->collSettings, &E);
	real_t D = evaluatePitchDistributionFromA(ir, xi0, A) * rGrid->GetVpVol(ir) / EvaluateVpREAtA(ir, A);
	return D;
}

/**
 * Evaluates the pitch distribution width parameter 'A'
 */
real_t AnalyticDistributionRE::GetAatP(len_t ir,real_t p, CollisionQuantity::collqty_settings *collSet, real_t *Eterm_in){
    CollisionQuantity::collqty_settings* settings = (collSet==nullptr) ? this->collSettings : collSet;
    real_t Eterm = (Eterm_in==nullptr) ? unknowns->GetUnknownData(id_Eterm)[ir] : *Eterm_in;
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(rGrid->GetFSA_B2(ir)); 
    real_t pNuD = p*nuD->evaluateAtP(ir,p,settings);    
    return 2*E/pNuD;
}

/**
 * Evaluate the full distribution function at the given electric field strength.
 */
real_t AnalyticDistributionRE::evaluateFullDistributionWithE(
	len_t ir, real_t xi0, real_t p, real_t E,
	real_t *dfdxi0, real_t *dfdp, real_t *dfdr
) {
	return evaluateEnergyDistribution(ir, p, dfdp, dfdr) *
		evaluatePitchDistributionWithE(ir, xi0, p, E, dfdxi0, dfdp, dfdr);
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
