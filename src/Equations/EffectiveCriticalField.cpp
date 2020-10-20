/**
 * Calculates the effective critical electric field
 * Settings parameters: 
 * - collQtySettings->collfreq_type (completely screened or not)
 * - Eceff_mode: 
 *  COLLQTY_ECEFF_MODE_EC_TOT = 0,      // Gives Ectot including all bound electrons (or Ec_free if no impurities/complete screening)
 *  COLLQTY_ECEFF_MODE_CYLINDRICAL = 1, // Sets Eceff using the Hesslow formula ignoring trapping effects.
 *  COLLQTY_ECEFF_MODE_SIMPLE = 2,      // An approximate numerical calculation with a simplified account of trapping effects
 *  COLLQTY_ECEFF_MODE_FULL = 3         // Full 'Lehtinen theory' expression.
 * 
 * 
 */
#include "FVM/config.h"
#include "DREAM/Equations/AnalyticDistribution.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/EffectiveCriticalField.hpp"


using namespace DREAM;

/**
 * Constructor.
 *
 */
EffectiveCriticalField::EffectiveCriticalField(ParametersForEceff *par)
    : Eceff_mode(par->Eceff_mode), collSettingsForEc(par->collSettingsForEc), collQtySettings(par->collQtySettings), 
    rGrid(par->rGrid), nuS(par->nuS), nuD(par->nuD),
    fsolve(par->fsolve)
    {
        gsl_parameters.rGrid = par->rGrid;
        gsl_parameters.nuS = par->nuS;
        gsl_parameters.nuD = par->nuD;
        gsl_parameters.fgType = par->fgType;
        gsl_parameters.gsl_ad_w = par->gsl_ad_w;
        gsl_parameters.fmin = par->fmin;
        gsl_parameters.collSettingsForEc = par->collSettingsForEc;
        gsl_parameters.QAG_KEY = GSL_INTEG_GAUSS31;
        gsl_parameters.analyticDist = new AnalyticDistribution(rGrid, nuD, Eceff_mode);
    }

EffectiveCriticalField::~EffectiveCriticalField(){
    delete gsl_parameters.analyticDist;
}


/**
 * Calculates and stores the effective critical field for runaway generation.
 * The calculation is based on Eq (21) in Hesslow et al, PPCF 60, 074010 (2018)
 * but has been generalized to account for inhomogeneous magnetic fields. The
 * method implemented here is outlined in DREAM/doc/notes/theory.
 * Essentially, the critical effective electric field is defined as the E
 * for which the maximum (with respect to p) of U(p) equals 0. Here, U
 * is the net momentum advection term averaged over an analytic pitch
 * angle distribution.
 */
void EffectiveCriticalField::CalculateEffectiveCriticalField(const real_t *Ec_tot, const real_t *Ec_free, real_t *effectiveCriticalField){
    len_t nr = rGrid->GetNr();
    // questions: 
    // why if-else-return here and not below (l.42 and onwards?) I think it makes more sense to have two returns or one big if statement 
    // should we account for radiation also if completely screened? it does with the other ECEFF_MODES, right?


    if(Eceff_mode == OptionConstants::COLLQTY_ECEFF_MODE_CYLINDRICAL){ 
        if(collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED)
            for(len_t ir=0; ir<nr; ir++)
                effectiveCriticalField[ir] = Ec_free[ir]; 
        else
            for(len_t ir=0; ir<nr; ir++)
                effectiveCriticalField[ir] = Ec_tot[ir];
        return;
    }
    // placeholder quantities that will be overwritten by the GSL functions
    std::function<real_t(real_t,real_t,real_t)> Func = [](real_t,real_t,real_t){return 0;};
    gsl_parameters.Func = Func; gsl_parameters.Eterm = 0; gsl_parameters.p = 0; gsl_parameters.p_ex_lo = 0;
    gsl_parameters.p_ex_up = 0;
    
    real_t ELo, EUp;
    gsl_function UExtremumFunc;
    for (len_t ir=0; ir<nr; ir++){
        gsl_parameters.ir = ir;
        UExtremumFunc.function = &(FindUExtremumAtE);
        UExtremumFunc.params = &gsl_parameters; // works with params here instead. 

        /**
         * Initial guess: Eceff is between 0.9*Ec_tot and 1.5*Ec_tot
         */
        ELo = .9*Ec_tot[ir];
        EUp = 1.5*Ec_tot[ir];
        RunawayFluid::FindInterval(&ELo, &EUp, UExtremumFunc);
        RunawayFluid::FindRoot(ELo,EUp, &effectiveCriticalField[ir], UExtremumFunc,fsolve);
    }
}

/**
 *  Returns the minimum of -U (with respect to p) at a given Eterm 
 */
real_t EffectiveCriticalField::FindUExtremumAtE(real_t Eterm, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    params->Eterm = Eterm;
    gsl_min_fminimizer *gsl_fmin = params->fmin;
    real_t p_ex_guess, p_ex_lo, p_ex_up;

    gsl_function F;
    F.function = &(UAtPFunc);
    F.params = params;
    real_t p_upper_threshold = 1000; // larger momenta are not physically relevant in our scenarios
    FindPExInterval(&p_ex_guess, &p_ex_lo, &p_ex_up, p_upper_threshold,params);

    // If the extremum is at a larger momentum than p_upper_threshold (or doesn't exist at all), 
    // we will define Eceff as the value where U(p_upper_threshold) = 0. 
    if(p_ex_up > p_upper_threshold)
        return UAtPFunc(p_upper_threshold,params);

    gsl_min_fminimizer_set(gsl_fmin, &F, p_ex_guess, p_ex_lo, p_ex_up);

    int status;
    real_t rel_error = 5e-2, abs_error=0;
    len_t max_iter = 30;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status     = gsl_min_fminimizer_iterate(gsl_fmin);
        p_ex_guess = gsl_min_fminimizer_x_minimum(gsl_fmin);
        p_ex_lo    = gsl_min_fminimizer_x_lower(gsl_fmin);
        p_ex_up    = gsl_min_fminimizer_x_upper(gsl_fmin);
        status     = gsl_root_test_interval(p_ex_lo, p_ex_up, abs_error, rel_error);

        if (status == GSL_SUCCESS)
            break;
    }

    real_t minimumFValue = gsl_min_fminimizer_f_minimum(gsl_fmin);
    return minimumFValue;
}


/**
 *  Finds an interval p \in [p_ex_lower, p_ex_upper] in which a minimum of -U(p) exists.  
 */
void EffectiveCriticalField::FindPExInterval(real_t *p_ex_guess, real_t *p_ex_lower, real_t *p_ex_upper, real_t p_upper_threshold,
                                            UContributionParams *params){
    *p_ex_lower = 1;
    *p_ex_upper = 100;
    *p_ex_guess = 10;
    real_t F_lo = UAtPFunc(*p_ex_lower,params);
    real_t F_up = UAtPFunc(*p_ex_upper,params);
    real_t F_g  = UAtPFunc(*p_ex_guess,params);
    
    if( (F_g < F_up) && (F_g < F_lo) ) // at least one minimum exists on the interval
        return;
    else if ( F_g > F_lo) // Minimum located at p<p_ex_guess
        while(F_g > F_lo){
            *p_ex_upper = *p_ex_guess;
            *p_ex_guess = *p_ex_lower;
            *p_ex_lower /= 5;
            F_g = F_lo; //UAtPFunc(*p_ex_guess,params);
            F_lo = UAtPFunc(*p_ex_lower,params);
        }
    else // Minimum at p>p_ex_guss
        while( (F_g > F_up) && (*p_ex_upper < p_upper_threshold)){
            *p_ex_lower = *p_ex_guess;
            *p_ex_guess = *p_ex_upper;
            *p_ex_upper *= 5;
            F_g = F_up;//UAtPFunc(*p_ex_guess,params);
            F_up = UAtPFunc(*p_ex_upper,params);
        }
}

/**
 * U(r,p) is the pitch-angle averaged momentum-advection coefficient,
 * U = < {A^p}f > / <f>,
 * where <...> = int ... d(xi0).
 * The calculation is documented under doc/notes/theory.pdf in 
 * Section 2 (under the heading 'Bounce-averaged effective field') 
 */

/*
The function takes a xi0 and a lambda expression Func (and other needed helper parameters) and 
returns the contribution to the integrand in the U function, i.e. V'{Func}*exp(-...),
where exp(-...)(xi0) is the analytical pitch-angle distribution, and V'{Func} the 
bounce integral of Func.
*/
real_t UPartialContribution(real_t xi0, void *par){
    struct EffectiveCriticalField::UContributionParams *params = (struct EffectiveCriticalField::UContributionParams *) par;
    CollisionQuantity::collqty_settings *collSettingsForEc = params->collSettingsForEc;
    FVM::RadialGrid *rGrid = params->rGrid; 
    len_t ir = params->ir;
    real_t p = params->p;
    FVM::fluxGridType fluxGridType = params->fgType;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    real_t E = params->Eterm;
    AnalyticDistribution *analyticDist = params-> analyticDist;
    std::function<real_t(real_t,real_t,real_t,real_t)> BAFunc = [xi0,params](real_t xiOverXi0,real_t BOverBmin,real_t /*ROverR0*/,real_t /*NablaR2*/){return params->Func(xi0,BOverBmin,xiOverXi0);};
    
    return rGrid->EvaluatePXiBounceIntegralAtP(ir,p,xi0,fluxGridType,BAFunc)
        * analyticDist->evaluatePitchDistribution(ir,xi0,p,E,collSettingsForEc, gsl_ad_w);
}

/**
 * Public method used mainly for benchmarking: evaluates the pitch-averaged friction function -U 
 */
real_t EffectiveCriticalField::testEvalU(len_t ir, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings){
    std::function<real_t(real_t,real_t,real_t)> Func = [](real_t,real_t,real_t){return 0;};
    real_t p_ex_lo = 0, p_ex_up = 0;
    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);
    const gsl_min_fminimizer_type *fmin_type = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *fmin = gsl_min_fminimizer_alloc(fmin_type);

    struct UContributionParams params = {gsl_parameters.rGrid, gsl_parameters.nuS, gsl_parameters.nuD, ir, p, 
    FVM::FLUXGRIDTYPE_DISTRIBUTION, Eterm, Func, gsl_ad_w,
                    fmin, p_ex_lo, p_ex_up,inSettings,gsl_parameters.QAG_KEY,gsl_parameters.analyticDist}; 
    return UAtPFunc(p,&params);
}

/**
 * Evaluates -U(p) at given Eterm.
 */
real_t EffectiveCriticalField::UAtPFunc(real_t p, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    params->p = p;
    FVM::RadialGrid *rGrid = params->rGrid;
    len_t ir = params->ir;
    FVM::fluxGridType fluxGridType = params->fgType;
    real_t Eterm = params->Eterm;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    SlowingDownFrequency *nuS = params->nuS;
    CollisionQuantity::collqty_settings *collSettingsForEc = params->collSettingsForEc;

    real_t Bmin,Bmax;
    if(fluxGridType == FVM::FLUXGRIDTYPE_RADIAL){
        Bmin = rGrid->GetBmin_f(ir);
        Bmax = rGrid->GetBmax_f(ir);    
    }else{
        Bmin = rGrid->GetBmin(ir);
        Bmax = rGrid->GetBmax(ir);
    }
    const real_t sqrtB2avgOverBavg = sqrt(rGrid->GetFSA_B2(ir)) / rGrid->GetFSA_B(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    if(xiT < 1e-6)
        xiT = 0;

    // Evaluates the contribution from electric field term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncElectric = 
            [](real_t xi0, real_t /*BOverBmin*/, real_t xiOverXi0 ){return xi0*xiOverXi0;};

    params->Func = FuncElectric;
    real_t EContrib, error;
    real_t Efactor = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrtB2avgOverBavg; 
    real_t epsabs = 0, epsrel = 5e-3, lim = gsl_ad_w->limit; 
    gsl_function GSL_func;
    GSL_func.function = &(UPartialContribution);
    GSL_func.params = params;
    if(xiT){
        real_t EContrib1, EContrib2;
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&EContrib1,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&EContrib2,&error);
        EContrib = EContrib1 + EContrib2;
    }else
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&EContrib,&error);
    EContrib *= Efactor;

    // Evaluates the contribution from slowing down term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncUnity = 
            [](real_t,real_t,real_t){return 1;};
    params->Func = FuncUnity;    
    real_t UnityContrib;
    if(xiT){
        real_t UnityContrib1, UnityContrib2, UnityContrib3;
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib1,&error);
        gsl_integration_qags(&GSL_func,-xiT,xiT,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib2,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib3,&error);
        UnityContrib = UnityContrib1 + UnityContrib2 + UnityContrib3;
    } else 
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib,&error);

    real_t NuSContrib = -p*nuS->evaluateAtP(ir,p,collSettingsForEc) * UnityContrib;

    // Evaluates the contribution from synchrotron term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncSynchrotron = 
            [](real_t xi0, real_t BOverBmin, real_t){return (1-xi0*xi0)*BOverBmin*BOverBmin*BOverBmin;};
    params->Func = FuncSynchrotron;
    real_t SynchrotronFactor = -p*sqrt(1+p*p)* Constants::ec * Constants::ec * Constants::ec * Constants::ec * Bmin * Bmin
                            / ( 6 * M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me
                                * Constants::c * Constants::c * Constants::c); 

    real_t SynchContrib;
    if(xiT){
        real_t SynchContrib1, SynchContrib2, SynchContrib3;
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib1,&error);
        gsl_integration_qags(&GSL_func,-xiT,xiT,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib2,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib3,&error);
        SynchContrib = SynchContrib1 + SynchContrib2 + SynchContrib3;
    } else 
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib,&error);

    SynchContrib *= SynchrotronFactor; 

    return -(EContrib + NuSContrib + SynchContrib) / UnityContrib;

}