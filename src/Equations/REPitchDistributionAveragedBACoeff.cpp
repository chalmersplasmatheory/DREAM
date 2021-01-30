
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
 */


REPitchDistributionAveragedBACoeff::REPitchDistributionAveragedBACoeff(
    FVM::RadialGrid *rg, AnalyticDistributionRE *dist,
    real_t(*func)(real_t,real_t,real_t,real_t,void*),
    void *par, int_t *list, 
    std::function<real_t(real_t)> pitchfunc,
    std::function<real_t(len_t,real_t)> momfunc,    
    const gsl_interp_type *t
) : rGrid(rg), distRE(dist), BA_Func(func), BA_Func_par(par), BA_Param(list), 
    BA_PitchPrefactor(pitchfunc), BA_MomentumPrefactor(momfunc), interp_mode(t) 
{
    gsl_ad_w = gsl_integration_workspace_alloc(1000);
    GridRebuilt();
}

REPitchDistributionAveragedBACoeff::~ REPitchDistributionAveragedBACoeff(){
    Deallocate();
    gsl_integration_workspace_free(gsl_ad_w);
}

bool REPitchDistributionAveragedBACoeff::GridRebuilt(){
    Deallocate();
    this->nr = rGrid->GetNr();

    generateBASplines();
    generateREDistAverageSplines();

    return true;
}


/**
 * Generates splines of BounceIntegral(BA_Func) on the interval
 * xi0 \in [0,1], since it by construction is a symmetric function
 * of xi0 (any asymmetries with respect to xi0 would be captured in
 * the 'PitchPrefactor' function instead)
 */
void REPitchDistributionAveragedBACoeff::generateBASplines(){
    real_t 
        xi0Array[N_BA_SPLINE],
        BAArray[N_BA_SPLINE];

    // xi0 uniform in [0,1]
    for(len_t i=0; i<N_BA_SPLINE; i++)
        xi0Array[i] = i*1.0/(N_BA_SPLINE-1);

    BA_Spline = new gsl_spline*[nr];
    BA_Accel  = new gsl_interp_accel*[nr];
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<N_BA_SPLINE; i++)
            BAArray[i]  = rGrid->EvaluatePXiBounceIntegralAtP(ir, xi0Array[i], FVM::FLUXGRIDTYPE_DISTRIBUTION, BA_Func, BA_Func_par, BA_Param);
        BA_Accel[ir] = gsl_interp_accel_alloc();
        BA_Spline[ir] = gsl_spline_alloc(interp_mode, N_BA_SPLINE);
        gsl_spline_init( BA_Spline[ir], xi0Array, BAArray, N_BA_SPLINE);
    }
}

// Integrand of the RE pitch distribution-bounce-integral.
struct ParametersForREPitchDistributionIntegral {len_t ir; real_t A; AnalyticDistributionRE *distRE; gsl_spline *spline; gsl_interp_accel *acc; std::function<real_t(real_t)> PitchFunc;};
real_t evaluateREPitchDistributionIntegralKernel(real_t xi0, void *par){
    ParametersForREPitchDistributionIntegral *params = (ParametersForREPitchDistributionIntegral*) par;
    len_t ir = params->ir;
    return params->PitchFunc(xi0)*gsl_spline_eval(params->spline, fabs(xi0), params->acc) 
        * params->distRE->evaluatePitchDistributionFromA(ir, xi0, params->A);
}

void REPitchDistributionAveragedBACoeff::generateREDistAverageSplines(){
    real_t error;
    real_t epsabs = 1e-6, epsrel = 5e-3, lim = gsl_ad_w->limit; 

    real_t 
        xArray[N_RE_DIST_SPLINE],
        REDistAverageArray[N_RE_DIST_SPLINE];

    for(len_t i=0; i<N_RE_DIST_SPLINE; i++)
        xArray[i] = i*1.0/(N_RE_DIST_SPLINE - 1);

    ParametersForREPitchDistributionIntegral params;
    params.distRE = distRE;
    gsl_function GSL_Func;
    GSL_Func.function = &(evaluateREPitchDistributionIntegralKernel);
    GSL_Func.params = &params;

    REDistAverage_Spline = new gsl_spline*[nr];
    REDistAverage_Accel  = new gsl_interp_accel*[nr];
    for(len_t ir=0; ir<nr; ir++){
        params.ir = ir;
        real_t xiT = rGrid->GetXi0TrappedBoundary(ir);

        params.spline = BA_Spline[ir];
        params.acc = BA_Accel[ir];
        params.PitchFunc = BA_PitchPrefactor;
        for(len_t i=0; i<N_RE_DIST_SPLINE; i++){
            real_t A = GetAFromX(xArray[i]);
            params.A = A;
            if(xiT){
                real_t contrib1, contrib2, contrib3;
                gsl_integration_qags(&GSL_Func, -1, -xiT, epsabs, epsrel, lim, gsl_ad_w, &contrib1, &error);
                gsl_integration_qags(&GSL_Func, 0, xiT, epsabs, epsrel, lim, gsl_ad_w, &contrib2, &error);
                gsl_integration_qags(&GSL_Func, xiT, 1, epsabs, epsrel, lim, gsl_ad_w, &contrib3, &error);
                REDistAverageArray[i] = (contrib1 + contrib2 + contrib3);
            } else 
                gsl_integration_qags(&GSL_Func, -1, 1, epsabs, epsrel, lim, gsl_ad_w, &REDistAverageArray[i], &error);

            REDistAverageArray[i] /= distRE->EvaluateVpREAtA(ir, A);
        }
        REDistAverage_Accel[ir] = gsl_interp_accel_alloc();
        REDistAverage_Spline[ir] = gsl_spline_alloc(interp_mode, N_RE_DIST_SPLINE);
        gsl_spline_init( REDistAverage_Spline[ir], xArray, REDistAverageArray, N_RE_DIST_SPLINE);
    }
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
