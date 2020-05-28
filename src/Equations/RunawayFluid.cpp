#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/NotImplementedException.hpp"
using namespace DREAM;

RunawayFluid::RunawayFluid(FVM::Grid *g, FVM::UnknownQuantityHandler *u, SlowingDownFrequency *nuS, 
    PitchScatterFrequency *nuD, CoulombLogarithm *lnLee, CollisionQuantity::collqty_settings *cqs){
    this->gridRebuilt = true;
    this->rGrid = g->GetRadialGrid();
    this->nuS = nuS;
    this->nuD = nuD;
    this->lnLambdaEE = lnLee;
    this->collQtySettings = cqs;
    this->unknowns = u;
    id_ncold = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ntot  = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_Tcold = this->unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Eterm = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);


    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    const gsl_min_fminimizer_type *fmin_type = gsl_min_fminimizer_brent;
    this->gsl_ad_w = gsl_integration_workspace_alloc(1000);
    this->fsolve = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    this->fmin = gsl_min_fminimizer_alloc(fmin_type);
}

RunawayFluid::~RunawayFluid(){
    DeallocateQuantities();

    gsl_integration_workspace_free(gsl_ad_w);
    gsl_root_fsolver_free(fsolve);
    gsl_min_fminimizer_free(fmin);
}

void RunawayFluid::Rebuild(bool useApproximateMethod){
    if(gridRebuilt){
        nr = rGrid->GetNr();
        AllocateQuantities();
    }
    ncold   = unknowns->GetUnknownData(id_ncold);
    ntot    = unknowns->GetUnknownData(id_ntot);
    Tcold   = unknowns->GetUnknownData(id_Tcold);
    Eterm   = unknowns->GetUnknownData(id_Eterm);
    
    CalculateDerivedQuantities();
    CalculateEffectiveCriticalField(useApproximateMethod);
    CalculateCriticalMomentum();
    CalculateGrowthRates();
}


void RunawayFluid::CalculateGrowthRates(){

    real_t *n_tot = unknowns->GetUnknownData(id_ntot); 
    real_t gamma_crit;
    for (len_t ir = 0; ir<this->nr; ir++){
        // we still haven't implemented the relativistic corrections in criticalREmomentum, 
        // but let's keep it like this for now in case we do in the future.

        if(criticalREMomentum[ir]==DBL_MAX){
            avalancheGrowthRate[ir] = 0;
        } else {
            gamma_crit = sqrt( 1 + criticalREMomentum[ir]*criticalREMomentum[ir] );
            avalancheGrowthRate[ir] = 0.5 * n_tot[ir] * constPreFactor / (gamma_crit-1) ;
        }
        tritiumRate[ir] = evaluateTritiumRate(gamma_crit);
        comptonRate[ir] = n_tot[ir]*evaluateComptonRate(criticalREMomentum[ir],gsl_ad_w);
    }
}


/**
 * Calculates the runaway rate due to beta decay of tritium. (to be multiplied by n_tritium)
 */
real_t RunawayFluid::evaluateTritiumRate(real_t gamma_c){
    real_t tau_halfLife = 12.32 * 365.24 *24*60*60; // 12.32 years, in seconds

    real_t decayMaxEnergyEV = 18.6e3; // maximum beta electron kinetic energy 
    real_t w = Constants::mc2inEV * (gamma_c-1) / decayMaxEnergyEV;
    real_t fracAbovePc = 1 + sqrt(w)*( -(35/8)*w + (21/4)*w*w - (15/8)*w*w*w);

    return log(2) /tau_halfLife * fracAbovePc;
}


// Evaluates total cross section for Compton scattering into p>pc due to incident photon of energy Eg (units of mc and mc2)
// Eq (29) in Martin-Solis NF 2017
real_t RunawayFluid::evaluateComptonTotalCrossSectionAtP(real_t Eg, real_t pc){
    real_t gamma_c = sqrt(1+pc*pc);
    real_t x = Eg;
    real_t Wc = pc*pc/(gamma_c+1); // = gamma_c-1
    real_t cc = 1 - 1/Eg * Wc /( Eg - Wc );
    return M_PI * Constants::r0 * Constants::r0 * ( (x*x-2*x-2)/(x*x*x) * log( (1+2*x)/( 1+x*(1-cc) ) ) 
        + 1/(2*x) * ( 1/( (1+x*(1-cc))*(1+x*(1-cc)) ) - 1/( (1+2*x)*(1+2*x) ) ) 
        - 1/(x*x*x) * ( 1 - x - (1+2*x) / (1+x*(1-cc)) - x*cc )   );
}

// Photon spectral flux density, Eq (24) in Martin-Solis NF 2017
real_t RunawayFluid::evaluateComptonPhotonFluxSpectrum(real_t Eg){
    real_t ITERPhotonFluxDensity = 1e18; // 1/m^2s
    real_t z = (1.2 + log(Eg * Constants::mc2inEV/1e6) ) / 0.8;
    return ITERPhotonFluxDensity * exp( - exp(z) - z + 1 );
}


// The integrand in the evaluation of the total production rate integral(flux density * cross section ) 
struct ComptonParam {real_t pc;};
real_t ComptonIntegrandFunc(real_t Eg, void *par){
    struct ComptonParam *params = (struct ComptonParam *) par;
    
    real_t pc = params->pc;

    return RunawayFluid::evaluateComptonPhotonFluxSpectrum(Eg) * RunawayFluid::evaluateComptonTotalCrossSectionAtP(Eg,pc);
}

// returns (dnRE/dt)_compton at radial index ir
real_t RunawayFluid::evaluateComptonRate(real_t pc,gsl_integration_workspace *gsl_ad_w){
    if(pc==DBL_MAX)
        return 0;

    real_t gamma_c = sqrt(1+pc*pc);
    real_t gammacMinusOne = pc*pc/(gamma_c+1); // = gamma_c-1
    struct ComptonParam  params= {pc};
    gsl_function ComptonFunc;
    ComptonFunc.function = &(ComptonIntegrandFunc);
    ComptonFunc.params = &params;

    real_t Eg_min = (pc + gammacMinusOne) /2;
    real_t valIntegral;
    // qagiu assumes an infinite upper boundary
    real_t epsrel = 1e-4;
    real_t epsabs;
    gsl_integration_qagiu(&ComptonFunc, Eg_min , 0, epsrel, 1000, gsl_ad_w, &valIntegral, &epsabs);
    return valIntegral;
}


struct pStarFuncParams {real_t constTerm; len_t ir; RunawayFluid *rf;};

// Returns gamma_trap^(1/4)*sqrt(E) * p - nuSbarnuDbar(p)^(1/4)
real_t RunawayFluid::pStarFunction(real_t p_eval, void *par){
    struct pStarFuncParams *params = (struct pStarFuncParams *) par;
    
    real_t constTerm = params->constTerm;
    real_t ir = params->ir;
    RunawayFluid *rf = params->rf;
    return sqrt(sqrt(rf->evaluateBarNuSNuDAtP(ir,p_eval))) - constTerm * p_eval ;

}
void RunawayFluid::CalculateCriticalMomentum(){

    real_t E, constTerm;
    real_t effectivePassingFraction;
    gsl_function gsl_func;
    pStarFuncParams pStar_params;
    real_t pLo, pUp, pStar;
    real_t *E_term = unknowns->GetUnknownData(id_Eterm); 
    for(len_t ir=0; ir<this->nr; ir++){
        if(E_term[ir] > effectiveCriticalField[ir])
            E =  Constants::ec * E_term[ir] /(Constants::me * Constants::c);
        else
            E =  Constants::ec * effectiveCriticalField[ir] /(Constants::me * Constants::c);

        
        /*
        Chooses whether trapping effects are accounted for in growth rates via setting 
        (could imagine another setting where you go smoothly from one to the other as 
        t_orbit/t_coll_at_pstar goes from <<1 to >>1)
        */
        if(collQtySettings->pstar_mode == OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONAL){
            effectivePassingFraction = 1;
        } else if(collQtySettings->pstar_mode == OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS){
            effectivePassingFraction = rGrid->GetEffPassFrac(ir);
        }
        constTerm = sqrt(sqrt(E*E * effectivePassingFraction));

        pStar_params = {constTerm,ir,this}; 
        
        gsl_func.function = &(pStarFunction);
        gsl_func.params = &pStar_params;

        real_t Ecfree_term = Constants::ec * Ec_free[ir] /(Constants::me * Constants::c);
        real_t Ectot_term  = Constants::ec * Ec_tot[ir]  /(Constants::me * Constants::c);
    
        pLo = 1/sqrt(sqrt(constTerm)/Ecfree_term);
        pUp = 1/sqrt(sqrt(pStar_params.constTerm)/Ectot_term);
        FindInterval(&pLo,&pUp, gsl_func);
        FindRoot(pLo,pUp, &pStar, gsl_func,fsolve);

        // Set critical RE momentum so that 1/critMom^2 = (E-Eceff)/sqrt(NuSbarNuDbar + 4*NuSbar)
        E = Constants::ec * (E_term[ir] - effectiveCriticalField[ir]) /(Constants::me * Constants::c);
        if (E<=0)
            criticalREMomentum[ir] = DBL_MAX; // should make growth rates zero
        else
            criticalREMomentum[ir] =  sqrt(sqrt( (evaluateBarNuSNuDAtP(ir,pStar) + 4*nuS->evaluateAtP(ir,pStar)*pStar*pStar*pStar/(1+pStar*pStar))  
                                                / (E*E * effectivePassingFraction) ));
    }
}



// Calculates the Connor-Hastie field Ec using the relativistic lnLambda
// and the Dreicer field ED using the thermal lnLambda.
void RunawayFluid::CalculateDerivedQuantities(){
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    for (len_t ir=0; ir<nr; ir++){
        Ec_free[ir] = lnLambdaEE->GetLnLambdaC(ir) * ncold[ir] * constPreFactor * Constants::me * Constants::c / Constants::ec;
        Ec_tot[ir]  = lnLambdaEE->GetLnLambdaC(ir) * ntot[ir]  * constPreFactor * Constants::me * Constants::c / Constants::ec;
        EDreic[ir]  = lnLambdaEE->GetLnLambdaT(ir) * ncold[ir] * constPreFactor * (Constants::me * Constants::c / Constants::ec) * (Constants::mc2inEV / T_cold[ir]);
        pc_COMPLETESCREENING[ir] = 1/sqrt(Eterm[ir]/Ec_free[ir] -1);
        pc_NOSCREENING[ir] = 1/sqrt(Eterm[ir]/Ec_tot[ir] -1);
    }
}



///////////////////////////////////////////////////////////////////
/////////// BEGINNING OF BLOCK WITH METHODS RELATED TO  ///////////
/////////// CALCULATION OF THE EFFECTIVE CRITICAL FIELD /////////// 
///////////////////////////////////////////////////////////////////
struct UContributionParams {FVM::RadialGrid *rGrid; /*FVM::MomentumGrid *mg;*/ RunawayFluid *rf; SlowingDownFrequency *nuS; PitchScatterFrequency *nuD; len_t ir; real_t p; bool rFluxGrid; 
                            real_t Eterm; std::function<real_t(real_t,real_t,real_t)> Func; gsl_integration_workspace *gsl_ad_w;
                            gsl_min_fminimizer *fmin;real_t p_ex_lo; real_t p_ex_up; bool useApproximateMethod;};



#include "RunawayFluid.UFuncApprox.cpp"
#include "RunawayFluid.UFuncExact.cpp"

real_t RunawayFluid::testEvalU(len_t ir, real_t p, real_t Eterm, bool useApproximateMethod){
    bool rFluxGrid = false;
    std::function<real_t(real_t,real_t,real_t)> Func = [](real_t,real_t,real_t){return 0;};
    real_t p_ex_lo = 0, p_ex_up = 0;
    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);
    const gsl_min_fminimizer_type *fmin_type = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *fmin = gsl_min_fminimizer_alloc(fmin_type);

    struct UContributionParams params = {rGrid, this, nuS,nuD, ir, p, rFluxGrid, Eterm, Func, gsl_ad_w,
                    fmin, p_ex_lo, p_ex_up,useApproximateMethod};
    return UAtPFunc(p,&params);
}


void RunawayFluid::CalculateEffectiveCriticalField(bool useApproximateMethod){
    effectiveCriticalField = new real_t[nr];
    
    // placeholder quantities that will be overwritten by the GSL functions
    std::function<real_t(real_t,real_t,real_t)> Func = [](real_t,real_t,real_t){return 0;};
    real_t Eterm = 0, p = 0, p_ex_lo = 0, p_ex_up = 0;

    bool rFluxGrid = false;
    real_t ELo, EUp;
    UContributionParams params; 
    gsl_function UExtremumFunc;
    for (len_t ir=0; ir<this->nr; ir++){
        params = {rGrid, this, nuS,nuD, ir, p, rFluxGrid, Eterm, Func, gsl_ad_w,
                            fmin, p_ex_lo, p_ex_up,useApproximateMethod};
        UExtremumFunc.function = &(FindUExtremumAtE);
        UExtremumFunc.params = &params;
        ELo = .9*Ec_tot[ir];
        EUp = 1.5*Ec_tot[ir];
        FindInterval(&ELo, &EUp, UExtremumFunc);
//        FindECritInterval(ir, &ELo, &EUp, &params);
        FindRoot(ELo,EUp, &effectiveCriticalField[ir], UExtremumFunc,fsolve);
    }
}


// Finds an E interval within which Eceff sits. It guesses that Eceff ~ Ectot and adjusts from there.
/*
void RunawayFluid::FindECritInterval(len_t ir, real_t *E_lower, real_t *E_upper, void *par){
     struct UContributionParams *params = (struct UContributionParams *) par;
     RunawayFluid *rf = params->rf;
    *E_lower = rf->Ec_tot[ir];

    // If E < Eceff, ECritFunc (the minimum of -U) will be positive, i.e. U=0 has no roots   
    bool isELoUnderestimate = (FindUExtremumAtE(*E_lower, params) > 0);
    while(!isELoUnderestimate){
        *E_upper = *E_lower;
        *E_lower *= 0.7;
        isELoUnderestimate = (FindUExtremumAtE(*E_lower, params) > 0);
    }
    if(!(*E_upper==0))
        return;

    *E_upper = 1.5*rf->Ec_tot[ir]; 
    bool isEUpOverestimate = (FindUExtremumAtE(*E_upper, params) < 0);
    while (!isEUpOverestimate){
        *E_lower = *E_upper;
        *E_upper *= 1.4;
        isEUpOverestimate = (FindUExtremumAtE(*E_upper, params) < 0);
    }

}
*/

// Returns the minimum of -U (with respect to p) at a given Eterm 
real_t RunawayFluid::FindUExtremumAtE(real_t Eterm, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    params->Eterm = Eterm;
    gsl_min_fminimizer *gsl_fmin = params->fmin;


    real_t p_ex_guess, p_ex_lo, p_ex_up;

    gsl_function F;
    F.function = &(UAtPFunc);
    F.params = params;
    real_t p_upper_threshold = 1000; // larger momenta are not physically relevant in our scenarios
    FindPExInterval(&p_ex_guess, &p_ex_lo, &p_ex_up, params, p_upper_threshold);

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

        if (status == GSL_SUCCESS){
            break;
        }
    }

    real_t minimumFValue = gsl_min_fminimizer_f_minimum(gsl_fmin);
    return minimumFValue;
}



// Finds an E interval within which Eceff sits. It guesses that Eceff ~ Ectot and adjusts from there.
void RunawayFluid::FindPExInterval(real_t *p_ex_guess, real_t *p_ex_lower, real_t *p_ex_upper, void *par, real_t p_upper_threshold){
    struct UContributionParams *params = (struct UContributionParams *) par;

    *p_ex_lower = 1;
    *p_ex_upper = 100;
    *p_ex_guess = 10;
    real_t F_lo = UAtPFunc(*p_ex_lower,params);
    real_t F_up = UAtPFunc(*p_ex_upper,params);
    real_t F_g  = UAtPFunc(*p_ex_guess,params);
    
    if( (F_g < F_up) && (F_g < F_lo) ) // at least one minimum exists on the interval
        return;
    else if ( F_g > F_lo){ // Minimum located at p<p_ex_guess
        while(F_g > F_lo){
            *p_ex_upper = *p_ex_guess;
            *p_ex_guess = *p_ex_lower;
            *p_ex_lower /= 5;
            F_g = F_lo; //UAtPFunc(*p_ex_guess,params);
            F_lo = UAtPFunc(*p_ex_lower,params);
        }
    } else { // Minimum at p>p_ex_guss
        while( (F_g > F_up) && (*p_ex_upper < p_upper_threshold)){
            *p_ex_lower = *p_ex_guess;
            *p_ex_guess = *p_ex_upper;
            *p_ex_upper *= 5;
            F_g = F_up;//UAtPFunc(*p_ex_guess,params);
            F_up = UAtPFunc(*p_ex_upper,params);
        }

    }

}

real_t RunawayFluid::UAtPFunc(real_t p, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    if (params->useApproximateMethod)
        return evaluateApproximateUAtP(p,par);
    else
        return evaluateNegUAtP(p,par);
}


// Takes a p interval [x_lower, x_upper] and iterates at most max_iter=30 times (or to a relative error of rel_error=0.001)
// to find an estimate for p_Star 
void RunawayFluid::FindRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func, gsl_root_fsolver *s){
    
    gsl_root_fsolver_set (s, &gsl_func, x_lower, x_upper); 

    int status;
    real_t epsrel = 3e-3;
    len_t max_iter = 30;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status   = gsl_root_fsolver_iterate (s);
        *root    = gsl_root_fsolver_root (s);
        x_lower = gsl_root_fsolver_x_lower (s);
        x_upper = gsl_root_fsolver_x_upper (s);
        status   = gsl_root_test_interval (x_lower, x_upper, 0, epsrel);

        if (status == GSL_SUCCESS){
            
            break;
        }
    }
}

// Sets the p interval for pStar to be between the completely screened and non-screened limits 
void RunawayFluid::FindInterval(real_t *x_lower, real_t *x_upper, gsl_function gsl_func ){
    // Guess: p_lower = completely screened pc
    //        p_upper = non-screened pc
    


    // If pStar is smaller than p_lower (for some reason), reduce by 30%
    bool isLoUnderestimate = (gsl_func.function(*x_lower, gsl_func.params ) > 0);
    bool isUpOverestimate = (gsl_func.function(*x_upper, gsl_func.params ) < 0);
    while(!isLoUnderestimate){
        *x_upper = *x_lower;
        *x_lower *= 0.7;
        isLoUnderestimate = (gsl_func.function(*x_lower, gsl_func.params ) > 0);
        isUpOverestimate = true;
    }
    while (!isUpOverestimate){
        *x_upper *= 1.4;
        isUpOverestimate = (gsl_func.function(*x_upper, gsl_func.params ) < 0);
    }
}






void RunawayFluid::AllocateQuantities(){
    DeallocateQuantities();
    Ec_free = new real_t[nr];
    Ec_tot  = new real_t[nr];
    EDreic  = new real_t[nr];
    effectiveCriticalField = new real_t[nr];
    criticalREMomentum     = new real_t[nr];
    pc_COMPLETESCREENING   = new real_t[nr];
    pc_NOSCREENING         = new real_t[nr];
    avalancheGrowthRate    = new real_t[nr];
    tritiumRate = new real_t[nr];
    comptonRate = new real_t[nr];



}

void RunawayFluid::DeallocateQuantities(){
    if(Ec_free != nullptr){
        delete [] Ec_free;
        delete [] Ec_tot;
        delete [] EDreic;
        delete [] effectiveCriticalField;
        delete [] criticalREMomentum;
        delete [] pc_COMPLETESCREENING;
        delete [] pc_NOSCREENING;
        delete [] avalancheGrowthRate;
        delete [] tritiumRate;
        delete [] comptonRate;

    }
}



