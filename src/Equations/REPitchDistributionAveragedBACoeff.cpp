
#include "DREAM/Equations/REPitchDistributionAveragedBACoeff.hpp"

using namespace DREAM;

/**
 * Implementation of a class which handles the averaging of
 * bounce averaged coefficients over the runaway pitch distribution,
 * in the form used in the calculation of the critical electric field.
 * 
 *      [[X]] = \int Vp {X} exp(-A*g(xi0)) dxi0 / VpRE = [X] / [1]
 *       VpRE = \int Vp exp(-A*g(xi0)) dxi0 \equiv [1]
 * 
 * here {X} is the regular bounce average, g the pitch-dependent exponent
 * appearing in the analytical runaway pitch distribution (from Lehtinen-style theory)
 * and [[X]] is the generalized bounce-pitch-average.
 * 
 * The instances of this class generate splines of the bounce integrals {X}*Vp/VpVol
 * as a function of pitch, as well as of the generalized bounce-distribution integrals [X]
 * as a function of the distribution width parameter 'A'.
 * 
 * The construction of the splines follows the following workflow:
 * 
 *  1) Construct a spline over the bounce integral of the unknown X
 *     as function of xi0 (on the reduced interval xi0 \in [0,1] since 
 *     the bounce integrals are always symmetric). 'GenerateBASpline(...)'.
 *  2) Integrate this spline weighted by the analytic distribution function
 *     and the provided pitch weight function 'BA_PitchPrefactor'. Normalize
 *     the result by the Vp-weighted pitch distribution integral and spline
 *     the result as a function of A (uniformly in a mapped variable X(A)).
 *     Carried out in 'generateREDistAverageSplines',
 *  3) The final averaged coefficient is obtained ('EvaluateREPitchDistAverage(...)') by 
 *     evaluating the spline from step 2), multiplied by the pitch-independent weight 
 *     function 'BA_MomentumPrefactor'.
 */


REPitchDistributionAveragedBACoeff::REPitchDistributionAveragedBACoeff(
    FVM::RadialGrid *rg, AnalyticDistributionRE *dist,
    real_t(*func)(real_t,real_t,real_t,real_t,void*),
    void *par, const int_t *list, 
    std::function<real_t(real_t)> pitchfunc,
    std::function<real_t(len_t,real_t)> momfunc,    
    const gsl_interp_type *t
) : rGrid(rg), distRE(dist), BA_Func(func), BA_Func_par(par), BA_Param(list), 
    BA_PitchPrefactor(pitchfunc), BA_MomentumPrefactor(momfunc), interp_mode(t) 
{
    gsl_ad_w = gsl_integration_workspace_alloc(1000);
}

REPitchDistributionAveragedBACoeff::~ REPitchDistributionAveragedBACoeff(){
    Deallocate();
    gsl_integration_workspace_free(gsl_ad_w);
}

bool REPitchDistributionAveragedBACoeff::GridRebuilt(){
    Deallocate();
    this->nr = rGrid->GetNr();

    BA_Spline = new gsl_spline*[nr];
    BA_Accel  = new gsl_interp_accel*[nr];
    for(len_t ir=0; ir<nr; ir++){
        BA_Accel[ir] = gsl_interp_accel_alloc();
        GenerateBASpline(
            ir, rGrid, rGrid->GetXi0TrappedBoundary(ir),
            N_BA_SPLINE, BA_Func, BA_Func_par, BA_Param,
            BA_Spline[ir], interp_mode
        );
    }
    generateREDistAverageSplines();

    return true;
}


/**
 * Creates a pitch grid which is spaced logarithmically around xiT,
 * uniform in |xi0-xiT| with points as close as xiT +/- exp(minArg).
 * In cylindrical plasmas the pitch grid is sampled uniformly
 */
void REPitchDistributionAveragedBACoeff::SetBASplineArray(real_t xiT, real_t *&xi0Array, len_t N, real_t fracPointsLower, real_t minArg){
    if(xiT==0){ // cylindrical plasma
        for(len_t i=0; i<N; i++)
            xi0Array[i] = i*1.0/(N-1.0);
        return;
    }

    len_t N1 = (len_t) (N * fracPointsLower); // rounds down to a natural number
    len_t N2 = N - N1;
    real_t maxArg = log(xiT);
    for(len_t i=0; i<N1; i++){
        real_t expFact = maxArg - i * (maxArg - minArg) / (N1-1.0); 
        xi0Array[i] = xiT - exp(expFact);
    }
    maxArg = log(1.0-xiT);
    for(len_t i=0; i<N2; i++){
        real_t expFact = maxArg - (N2-1.0-i) * (maxArg - minArg) / (N2-1.0); 
        xi0Array[N1 + i] = xiT + exp(expFact);
    }
}


/**
 * Generates splines of BounceIntegral(BA_Func) on the interval
 * xi0 \in [0,1], since it by construction is a symmetric function
 * of xi0 (any asymmetries with respect to xi0 would be captured in
 * the 'PitchPrefactor' function instead)
 */
void REPitchDistributionAveragedBACoeff::GenerateBASpline(
    len_t ir, FVM::RadialGrid *rGrid, real_t xiT, const len_t N, 
    real_t(*Func)(real_t,real_t,real_t,real_t,void*),
    void *Func_par, const int_t *Param, 
    gsl_spline *&spline, const gsl_interp_type *interp_mode, 
    real_t fracPointsLower, real_t minArg
){
    real_t 
        *xi0Array = new real_t[N],
        *BAArray  = new real_t[N];

    SetBASplineArray(xiT, xi0Array, N, fracPointsLower, minArg);
    for(len_t i=0; i<N; i++)
        BAArray[i]  = rGrid->EvaluatePXiBounceIntegralAtP(
            ir, xi0Array[i], FVM::FLUXGRIDTYPE_DISTRIBUTION, 
            Func, Func_par, Param
        );
    spline = gsl_spline_alloc(interp_mode, N);
    gsl_spline_init( spline, xi0Array, BAArray, N);

    delete [] xi0Array;
    delete [] BAArray;
}


/**
 * Create xArray on [0,1], on which we will interpolate 
 * the RE pitch distribution averaged coefficients.
 * Constructed as:
 *  N * fracPointsLower on the interval [0,1-fracUpperInterval]
 *  remaining points on the interval [1-fracUpperInterval, 1]
 * corresponding (here) to a denser grid near xArray ~ 1, corresponding to
 * large A (strong electric fields, beam-like distributions)
 */
void REPitchDistributionAveragedBACoeff::GenerateNonUniformXArray(
    real_t *&xArray, const len_t N, real_t fracPointsLower, real_t fracUpperInterval
){
    len_t N1 = (len_t) (N * fracPointsLower); // rounds down to a natural number
    len_t N2 = N - N1;
    for(len_t i=0; i<N1; i++)
        xArray[i] = i*(1.0-fracUpperInterval)/(N1-1.0);
    for(len_t i=1; i<=N2; i++)
        xArray[N1-1+i] = 1.0 - fracUpperInterval + i*fracUpperInterval/N2;
}

void REPitchDistributionAveragedBACoeff::generateREDistAverageSplines(){
    real_t 
        *xArray = new real_t[N_RE_DIST_SPLINE],
        *REDistAverageArray = new real_t[N_RE_DIST_SPLINE];

    GenerateNonUniformXArray(xArray, N_RE_DIST_SPLINE);

    REDistAverage_Spline = new gsl_spline*[nr];
    REDistAverage_Accel  = new gsl_interp_accel*[nr];
    for(len_t ir=0; ir<nr; ir++){
        REDistAverage_Accel[ir] = gsl_interp_accel_alloc();
        REDistAverage_Spline[ir] = gsl_spline_alloc(interp_mode, N_RE_DIST_SPLINE);
        real_t xiT = rGrid->GetXi0TrappedBoundary(ir);
        ParametersForREPitchDistributionIntegral params = 
            {ir, xiT, 0, distRE, BA_Spline[ir], BA_Accel[ir], BA_PitchPrefactor};
        for(len_t i=0; i<N_RE_DIST_SPLINE-1; i++){
            real_t A = GetAFromX(xArray[i]);
            params.A = A;
            REDistAverageArray[i] = EvaluateREDistBounceIntegral(params, gsl_ad_w) 
                                    / distRE->EvaluateVpREAtA(ir, A);
        }
        // following two calls: set the singular point A=inf with a reduced form applicable to this limit
        real_t BAAtUnityXi = rGrid->CalculatePXiBounceAverageAtP(
            ir, 1.0, FVM::FLUXGRIDTYPE_DISTRIBUTION, 
            BA_Func, BA_Func_par, BA_Param
        );
        REDistAverageArray[N_RE_DIST_SPLINE-1] = params.PitchFunc(1.0)*BAAtUnityXi; //gsl_spline_eval(params.spline, 1.0, params.acc);

        gsl_spline_init(REDistAverage_Spline[ir], xArray, REDistAverageArray, N_RE_DIST_SPLINE);
    }
    delete [] xArray;
    delete [] REDistAverageArray;
}

/**
 * Integrand of the RE pitch distribution-bounce-integral.
 */
real_t REPitchDistributionAveragedBACoeff::evaluateREPitchDistributionIntegralKernel(real_t xi0, void *par){
    ParametersForREPitchDistributionIntegral *params = (ParametersForREPitchDistributionIntegral*) par;
    if(xi0<=0 && xi0>=-params->xiT)
        return 0;
    return params->PitchFunc(xi0)*gsl_spline_eval(params->spline, fabs(xi0), params->acc) 
        * params->distRE->evaluatePitchDistributionFromA(params->ir, xi0, params->A);
}

/**
 * Integrates the non-normalized pitch distribution over pitch, weighted 
 * with the params.PitchFunc function.
 */
real_t REPitchDistributionAveragedBACoeff::EvaluateREDistBounceIntegral(
    ParametersForREPitchDistributionIntegral params, gsl_integration_workspace *gsl_ad_w
){
    /**
     * Assumes that the non-normalized pitch distribution integrates 
     * to 0 in the limit of an infinitely narrow distribution
     */
    if(isinf(params.A))
        return 0.0;
    
    real_t xiT = params.xiT;

    real_t error;
    real_t epsabs = 1e-6, epsrel = 1e-3, lim = gsl_ad_w->limit; 
    static constexpr len_t numQAGPBreakpoints = 5; // integrate from -1 to 1 with QAGP
    real_t QAGPBreakpoints[numQAGPBreakpoints] = {-1.0, -xiT, 0.0, xiT, 1.0};

    gsl_function GSL_Func;
    GSL_Func.function = &(evaluateREPitchDistributionIntegralKernel);
    GSL_Func.params = &params;
    
    real_t val;
    if(xiT)
        gsl_integration_qagp(&GSL_Func, QAGPBreakpoints, numQAGPBreakpoints, epsabs, epsrel, lim, gsl_ad_w, &val, &error);
    else 
        gsl_integration_qags(&GSL_Func, -1.0, 1.0, epsabs, epsrel, lim, gsl_ad_w, &val, &error);

    return val;
}

/**
 * Main function of this class: evaluates 
 * the distribution-bounce averaged function, [[X]] 
 * at pitch distribution width parameter 'A'.
 * For optimization purposes, A can be provided
 * as an input, but the momentum 'p' will not be 
 * adjusted to be consistent with the provided value.
 */
real_t REPitchDistributionAveragedBACoeff::EvaluateREPitchDistAverage(len_t ir, real_t p, real_t *A_in){
    real_t preFactor = BA_MomentumPrefactor(ir,p);
    real_t A = (A_in == nullptr) ? distRE->GetAatP(ir,p) : *A_in;
    real_t X = GetXFromA(A);
    real_t distAverage = gsl_spline_eval(REDistAverage_Spline[ir], X, REDistAverage_Accel[ir]);
    real_t Y = preFactor * distAverage;
    
    /*if(dYdp != nullptr){ TODO: p derivative (if need be)
        dYdp = preFactor * dAdp * dXdA * gsl_spline_eval_deriv(REDistAverage_Spline[ir], X, REDistAverage_Accel[ir])
            + DDPpreFactor * distAverage;
        
        throw DREAMException(
            "REPitchDistributionAveragedCoeff: Evaluation of the p derivative" 
            "of the averaged coefficients has not yet been implemented.");
    }*/
    return Y;
}

/**
 * Deallocator
 */
void REPitchDistributionAveragedBACoeff::Deallocate(){
    if(BA_Spline != nullptr){
        for(len_t ir=0; ir<nr; ir++){
            gsl_spline_free(BA_Spline[ir]);
            gsl_interp_accel_free(BA_Accel[ir]);
        }
        delete [] BA_Spline;
        delete [] BA_Accel;
    }
    if(REDistAverage_Spline != nullptr){
        for(len_t ir=0; ir<nr; ir++){
            gsl_spline_free(REDistAverage_Spline[ir]);
            gsl_interp_accel_free(REDistAverage_Accel[ir]);
        }
        delete [] REDistAverage_Spline;
        delete [] REDistAverage_Accel;
    }
}
