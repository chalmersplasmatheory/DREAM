/**
 * Implementation of collision-rate calculator that calculates
 * various collision rates and related quantities, such as runaway growth rates.
*/



#include "DREAM/Equations/CollisionQuantityHandler.hpp"
//#include "DREAM/Constants.hpp"
//#include "DREAM/Settings/OptionConstants.hpp"
//#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
//#include <iostream>
//#include <fstream>


using namespace DREAM;

const len_t  CollisionQuantityHandler::conductivityLenT = 14;
const len_t  CollisionQuantityHandler::conductivityLenZ = 6;
const real_t CollisionQuantityHandler::conductivityBraams[conductivityLenZ*conductivityLenT] = {3.75994, 3.7549, 3.7492, 3.72852, 3.6842, 3.57129, 3.18206, 2.65006, 2.03127, 1.33009, 0.94648, 0.67042, 0.42422, 0.29999, 7.42898, 7.27359, 7.12772, 6.73805, 6.20946, 5.43667, 4.13733, 3.13472, 2.27862, 1.45375, 1.02875, 0.72743, 0.46003, 0.32528, 8.7546, 8.53281, 8.32655, 7.78445, 7.06892, 6.06243, 4.47244, 3.32611, 2.39205, 1.51805, 1.07308, 0.75853, 0.47965, 0.33915, 10.39122, 10.07781, 9.78962, 9.04621, 8.09361, 6.80431, 4.8805, 3.57303, 2.54842, 1.61157, 1.13856, 0.80472, 0.50885, 0.35979, 11.33006, 10.95869, 10.61952, 9.75405, 8.66306, 7.21564, 5.11377, 3.72206, 2.64827, 1.67382, 1.18263, 0.83593, 0.52861, 0.37377, 12.76615, 12.29716, 11.87371, 10.81201, 9.50746, 7.82693, 5.47602, 3.96944, 2.82473, 1.7887, 1.2649, 0.89443, 0.56569, 0.4};
const real_t CollisionQuantityHandler::conductivityTmc2[conductivityLenT] = {0,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100};
const real_t CollisionQuantityHandler::conductivityX[conductivityLenZ]    = {0,0.090909090909091,0.166666666666667,0.333333333333333,0.5,1};

/** 
 * Constructor
 */ 
CollisionQuantityHandler::CollisionQuantityHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantity::collqty_settings *cqset){
    ionHandler = ih;
    grid       = g;
    unknowns   = u;
    settings   = cqset;
    gridtype   = mgtype;

    this->nr   = g->GetNr();

    lnLambdaEE = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, settings,CollisionQuantity::LNLAMBDATYPE_EE);
    lnLambdaEI = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, settings,CollisionQuantity::LNLAMBDATYPE_EI);
    nuS   = new SlowingDownFrequency(grid, unknowns, ionHandler, lnLambdaEE,lnLambdaEI,gridtype, settings);
    nuD   = new PitchScatterFrequency(grid, unknowns, ionHandler, lnLambdaEI,lnLambdaEE,gridtype, settings);
    nuPar = new ParallelDiffusionFrequency(grid, unknowns, ionHandler, nuS,lnLambdaEE, gridtype, settings);
    REFluid = new RunawayFluid(grid,unknowns,nuS,nuD,lnLambdaEE);

    gsl_interp2d_init(gsl_cond, conductivityTmc2, conductivityX, conductivityBraams,conductivityLenT,conductivityLenZ);
}

/**
 * Destructor.
 */
CollisionQuantityHandler::~CollisionQuantityHandler(){
    DeallocateDerivedQuantities();
    delete [] nuS;
    delete [] nuD;
    delete [] nuPar;
    delete [] lnLambdaEE;
    delete [] lnLambdaEI;

    gsl_interp2d_free(gsl_cond);
}

/**
 * Calculates and stores all collision quantities from an UnknownQuantityHandler. 
 */
void CollisionQuantityHandler::Rebuild() {

    len_t id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->n_cold   = unknowns->GetUnknownData(id_ncold);

    len_t id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->T_cold   = unknowns->GetUnknownData(id_Tcold);

    len_t id_ntot  = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    this->n_tot    = unknowns->GetUnknownData(id_ntot);

    len_t id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->E_term   = unknowns->GetUnknownData(id_Eterm);

    this->nZ = ionHandler->GetNZ();
    this->nzs = ionHandler->GetNzs();
    this->ZAtomicNumber = ionHandler->GetZs();

    RunawayFluid *REFluid2 = new RunawayFluid(grid,unknowns,nuS,nuD,lnLambdaEE);
    lnLambdaEE->Rebuild();
    lnLambdaEI->Rebuild();
    nuS->Rebuild();
    nuD->Rebuild();
    nuPar->Rebuild();
    REFluid2->Rebuild(false);
    REFluid->Rebuild(true);
    CalculateDerivedQuantities();

}

void CollisionQuantityHandler::gridRebuilt(){
    lnLambdaEE->GridRebuilt();
    lnLambdaEI->GridRebuilt();
    nuS->GridRebuilt();
    nuD->GridRebuilt();
    nuPar->GridRebuilt();
}

// Uses collision frequencies and ion species to calculate
// critical fields and avalanche growth rates
void CollisionQuantityHandler::CalculateDerivedQuantities(){
    DeallocateDerivedQuantities();

    Ec_free = new real_t[this->nr];
    Ec_tot  = new real_t[this->nr];
    EDreic  = new real_t[this->nr];

    // Calculates the Connor-Hastie field Ec using the relativistic lnLambda
    // and the Dreicer field ED using the thermal lnLambda.
    for (len_t ir=0; ir<this->nr; ir++){
        Ec_free[ir] = lnLambdaEE->GetLnLambdaC(ir) * n_cold[ir] * constPreFactor * Constants::me * Constants::c / Constants::ec;
        Ec_tot[ir]  = lnLambdaEE->GetLnLambdaC(ir) * n_tot[ir]  * constPreFactor * Constants::me * Constants::c / Constants::ec;
        EDreic[ir]  = lnLambdaEE->GetLnLambdaT(ir) * n_cold[ir] * constPreFactor * (Constants::me * Constants::c / Constants::ec) * (Constants::mc2inEV / T_cold[ir]);
    }
    CalculateEffectiveCriticalField();
    CalculatePStar();

    CalculateGrowthRates();
}

real_t CollisionQuantityHandler::evaluateElectricalConductivity(len_t ir){
    const real_t T_SI = T_cold[ir] * Constants::ec;
    const real_t *Zeff = ionHandler->evaluateZeff();

    real_t sigmaBar = gsl_interp2d_eval(gsl_cond, conductivityTmc2, conductivityX, conductivityBraams, 
                T_SI / (Constants::me * Constants::c * Constants::c), 1/(1+Zeff[ir]), gsl_xacc, gsl_yacc  );
    
    real_t BraamsConductivity = 4*M_PI*Constants::eps0*Constants::eps0 * T_SI*sqrt(T_SI) / 
            (sqrt(Constants::me) * Constants::ec * Constants::ec * lnLambdaEE->GetLnLambdaT(ir) ) * sigmaBar;
    delete [] Zeff;
    return BraamsConductivity;
}



void CollisionQuantityHandler::DeallocateDerivedQuantities(){
    if (this->Ec_free == nullptr)
        return;

    delete [] this->Ec_free;
    delete [] this->Ec_tot;
    delete [] this->EDreic;
    delete [] this->criticalREMomentum;
    delete [] this->avalancheRate;
    delete [] this->tritiumRate;
    delete [] this->comptonRate;
    delete [] this->effectiveCriticalField;
}




// For GSL functions: partial contributions to evaluateUAtP
struct UFuncParams {len_t ir; real_t A; FVM::RadialGrid *rGrid;};
real_t UFrictionTermIntegrandCQH(real_t xi0, void *par){
    struct UFuncParams *params = (struct UFuncParams *) par;
    len_t ir = params->ir;
    real_t A = params->A;
    FVM::RadialGrid *rGrid = params->rGrid;
    if(xi0==0)
        return exp(-A);
    std::function<real_t(real_t,real_t,real_t)> FrictionTermFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return xi0/sqrt(1-BOverBmin*(1-xi0*xi0))*BOverBmin;};
    return rGrid->CalculateFluxSurfaceAverage(ir,false, FrictionTermFunc)*exp(-A*(1-xi0));
}

real_t USynchrotronTermIntegrandCQH(real_t xi0, void *par){
    struct UFuncParams *params = (struct UFuncParams *) par;
    len_t ir = params->ir;
    real_t A = params->A;
    FVM::RadialGrid *rGrid = params->rGrid;
    if(xi0==0)
        return exp(-A);
    std::function<real_t(real_t,real_t,real_t)> SynchrotronTermFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return (1-xi0*xi0)*xi0/sqrt(1-BOverBmin*(1-xi0*xi0)) *BOverBmin*BOverBmin*BOverBmin*BOverBmin;};
    return rGrid->CalculateFluxSurfaceAverage(ir,false, SynchrotronTermFunc)*exp(-A*(1-xi0));
}


// Evaluates the effective momentum flow U accounting for electric field, collisional friction and radiation reaction 
real_t CollisionQuantityHandler::evaluateUAtP(len_t ir,real_t p, real_t Eterm,gsl_integration_workspace *gsl_ad_w){
    FVM::RadialGrid *rGrid =  grid->GetRadialGrid();
    const real_t Bmin = rGrid->GetBmin(ir);
    const real_t Bmax = rGrid->GetBmax(ir);
    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 
    real_t xiT = sqrt(1-Bmin/Bmax);
//    real_t xiT = -1;
    real_t pNuD = p*nuD->evaluateAtP(ir,p,settings->collfreq_type, OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL);
    real_t A = 2*E/pNuD;
    real_t Econtrib;
    if(A==0)
        Econtrib = 0.5*E*(1-xiT*xiT);
    else 
        Econtrib = E/(A*A) *( A-1 - exp(-A*(1-xiT))*(A*xiT -1) );

    real_t FrictionTerm = p*nuS->evaluateAtP(ir,p,settings->collfreq_type, OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL);
    if(!(settings->bremsstrahlung_mode==OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT))
        FrictionTerm += evaluateBremsStoppingForceAtP(ir,p) / (Constants::me * Constants::c);
    
    UFuncParams FuncParams = {ir, A, rGrid};
    gsl_function UIntegrandFunc;

    UIntegrandFunc.function = &(UFrictionTermIntegrandCQH);
    UIntegrandFunc.params = &FuncParams;
    real_t abserr;
    real_t frictionIntegral;
//    real_t frictionIntegralPredict = (1-exp(-A))/A;
    gsl_integration_qags(&UIntegrandFunc, xiT,1.0,0,1e-4,1000,gsl_ad_w, &frictionIntegral, &abserr);
    real_t FrictionContrib = -FrictionTerm * frictionIntegral;

    real_t SynchrotronTerm = p*sqrt(1+p*p)* Constants::ec * Constants::ec * Constants::ec * Constants::ec * Bmin * Bmin
                            / ( 6 * M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me
                                * Constants::c * Constants::c * Constants::c); 
                               
    UIntegrandFunc.function = &(USynchrotronTermIntegrandCQH);
    real_t synchrotronIntegral;
//    real_t synchrotronIntegralPredict = exp(-A)*(-A*A+2*(A-1)*exp(A)+2)/(A*A*A);
    gsl_integration_qags(&UIntegrandFunc, xiT,1.0,0,1e-4,1000,gsl_ad_w, &synchrotronIntegral, &abserr);
    real_t SynchrotronContrib = -SynchrotronTerm * synchrotronIntegral;
    
   return (Econtrib + FrictionContrib + SynchrotronContrib) / frictionIntegral;
}

// For GSL function: returns U(p) at a given Eterm -- This could be compactified by writing evaluateUAtP on this form from the beginning
real_t UExtremumFunc(real_t p, void *par){
    struct CollisionQuantityHandler::UExtremumParams *params = (struct CollisionQuantityHandler::UExtremumParams *) par;
    len_t ir = params->ir;
    real_t Eterm = params->Eterm;
    gsl_integration_workspace *gsl_w = params->gsl_w;
    CollisionQuantityHandler *collQtyHand = params->collQtyHand;
    return - collQtyHand->evaluateUAtP(ir,p,Eterm,gsl_w);
}

// Returns the minimum of -U (with respect to p) at a given Eterm 
real_t FindUExtremumAtE(len_t ir, real_t Eterm, real_t *p_ex, gsl_integration_workspace *gsl_ad_w,CollisionQuantityHandler *collQtyHand){
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *gsl_fmin;

    // Requiring that the solution lies between p=0.1 and p=1000... 
    // should write a function that estimates these three based on non-screened and completely screened limits or something
    real_t p_ex_guess = 10.0;
    real_t p_ex_lo = 0.001, p_ex_up = 1000.0;
    gsl_function F;

    CollisionQuantityHandler::UExtremumParams params = {ir,Eterm,gsl_ad_w,collQtyHand};
    F.function = &(UExtremumFunc);
    F.params = &params;

    // TODO: Move allocation of s to CalculateEffectiveCriticalField()    
    T = gsl_min_fminimizer_brent;
    gsl_fmin = gsl_min_fminimizer_alloc(T);
    gsl_min_fminimizer_set(gsl_fmin, &F, p_ex_guess, p_ex_lo, p_ex_up);


    int status;
    real_t rel_error = 1e-2;
    len_t max_iter = 30;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status  = gsl_min_fminimizer_iterate(gsl_fmin);
        *p_ex   = gsl_min_fminimizer_x_minimum(gsl_fmin);
        p_ex_lo = gsl_min_fminimizer_x_lower(gsl_fmin);
        p_ex_up = gsl_min_fminimizer_x_upper(gsl_fmin);
        status  = gsl_root_test_interval(p_ex_lo, p_ex_up, 0, rel_error);

        if (status == GSL_SUCCESS){
            break;
        }
    }
/*
    std::cout << "Eterm:" << Eterm << std::endl;
    std::cout << "p_min:" << *p_ex << std::endl << std::endl;
*/  
    real_t minimumFValue = gsl_min_fminimizer_f_minimum(gsl_fmin);
    gsl_min_fminimizer_free(gsl_fmin);
    // Return function value -U(p_ex)
    return minimumFValue;
}

// For GSL function: returns minimum of -U at given Eterm -- This could be compactified by writing FindUExtremumAtE on this form from the beginning
real_t ECritFunc(real_t E, void *par){
    struct CollisionQuantityHandler::UExtremumParams *params = (struct CollisionQuantityHandler::UExtremumParams *) par;
    len_t ir = params->ir;
    gsl_integration_workspace *gsl_w = params->gsl_w;
    CollisionQuantityHandler *collQtyHand = params->collQtyHand;

    real_t p_extremum;
    return FindUExtremumAtE(ir, E, &p_extremum, gsl_w,collQtyHand);    
}


/** 
 * The function evaluates the effective critical field, defined (in doc/notes/theory.pdf) as the minimum electric field
 * (or rather, Eterm, effectively a loop voltage) above which U(p)=0 has real roots in p, where U is the pitch
 * averaged momentum flux in the limit (E-Eceff)<<E. It can be a bit hard to penetrate due to all GSL function thrown around,
 * but essentially it performs a root finding algorithm to solve the problem U_max(Eceff) = 0. U_max is in turn obtained using
 * a minimization algorithm in -U(p; E), where finally U(p) is obtained using gsl integration of various flux surface averages.
 * It performs up to 900 (max_iter*max_iter) function evaluations of nu_s and nu_D for each radial index ir, and in each such evaluation also evaluates 
 * multiple flux surface averages (inside a gsl adaptive quadrature, so tens?).
 */
void CollisionQuantityHandler::CalculateEffectiveCriticalField(){
    effectiveCriticalField = new real_t[nr];
    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);
    
    //real_t Eceff_guess;

    real_t ELo, EUp;
    UExtremumParams params;
    gsl_function UExtremumFunc;
    for (len_t ir=0; ir<this->nr; ir++){
        params = {ir,0,gsl_ad_w,this};
        UExtremumFunc.function = &(ECritFunc);
        UExtremumFunc.params = &params;

        ELo = 0;
        EUp = 0;
        FindECritInterval(ir, &ELo, &EUp, params);
        FindPStarRoot(ELo,EUp, &effectiveCriticalField[ir], UExtremumFunc);
    }

    gsl_integration_workspace_free(gsl_ad_w);

    /**
     * 1) Use optimization algorithm to find p_ex which minimizes -U(p; Eterm) at given Eterm
     * 2) Use root finding algorithm to find the Eterm for which U(p_ex; Eterm) = 0
     * 3) ???
     * 4) profit.
     */
}

// Finds an E interval within which Eceff sits. It guesses that Eceff ~ Ectot and adjusts from there.
void CollisionQuantityHandler::FindECritInterval(len_t ir, real_t *E_lower, real_t *E_upper, UExtremumParams params){

    *E_lower = Ec_tot[ir];

    // If E < Eceff, ECritFunc (the minimum of -U) will be positive, i.e. U=0 has no roots   
    bool isELoUnderestimate = (ECritFunc(*E_lower, &params) > 0);
    while(!isELoUnderestimate){
        *E_upper = *E_lower;
        *E_lower *= 0.7;
        isELoUnderestimate = (ECritFunc(*E_lower, &params) > 0);
    }
    if(!(*E_upper==0))
        return;

    *E_upper = 1.5*Ec_tot[ir]; 
    bool isEUpOverestimate = (ECritFunc(*E_upper, &params) < 0);
    while (!isEUpOverestimate){
        *E_lower = *E_upper;
        *E_upper *= 1.4;
        isEUpOverestimate = (ECritFunc(*E_upper, &params) < 0);
    }

}




// Returns gamma_trap^(1/4)*sqrt(E) * p - nuSbarnuDbar(p)^(1/4)
real_t CollisionQuantityHandler::pStarFunction(real_t p_eval, void *par){
    struct pStarFuncParams *params = (struct pStarFuncParams *) par;
    
    real_t constTerm = params->constTerm;
    real_t ir = params->ir;
    CollisionQuantityHandler *collQtyHand = params->collQtyHand;
    return constTerm * p_eval - sqrt(sqrt(collQtyHand->evaluateBarNuSNuDAtP(ir,p_eval)));

}

/**
 * Evaluates criticalREMomentum based on pStar which satisfies
 * pStar^2*nu_s(pStar)*nu_D(pStar) = ec*ec* Eterm*Eterm * effectivePassingFraction.
 */
void CollisionQuantityHandler::CalculatePStar(){
    criticalREMomentum = new real_t[this->nr];


    real_t E, constTerm;
    real_t effectivePassingFraction;
    gsl_function gsl_func;
    pStarFuncParams pStar_params;
    real_t pLo, pUp, pStar;
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
        if(settings->pstar_mode == OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONAL){
            effectivePassingFraction = 1;
        } else if(settings->pstar_mode == OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS){
            effectivePassingFraction = grid->GetRadialGrid()->GetEffPassFrac(ir);
        }
        constTerm = sqrt(sqrt(E*E * effectivePassingFraction));

        pStar_params = {constTerm,ir,this}; 
        
        gsl_func.function = &(pStarFunction);
        gsl_func.params = &pStar_params;


        pLo = 0;
        pUp = 0;
        FindPInterval(ir,&pLo,&pUp, pStar_params);
        pStar = 0;
        FindPStarRoot(pLo,pUp, &pStar, gsl_func);

        // Set critical RE momentum so that 1/critMom^2 = (E-Eceff)/sqrt(NuSbarNuDbar + 4*NuSbar)
        E = Constants::ec * (E_term[ir] - effectiveCriticalField[ir]) /(Constants::me * Constants::c);
        if (E<=0)
            criticalREMomentum[ir] = DBL_MAX;
        else
            criticalREMomentum[ir] =  sqrt(sqrt( (evaluateBarNuSNuDAtP(ir,pStar) + 4*nuS->evaluateAtP(ir,pStar)*pStar*pStar*pStar/(1+pStar*pStar))  
                                                / (E*E * effectivePassingFraction) ));
    }
}

// Sets the p interval for pStar to be between the completely screened and non-screened limits 
void CollisionQuantityHandler::FindPInterval(len_t ir, real_t *p_lower, real_t *p_upper, pStarFuncParams pStar_params ){
    // Guess: p_lower = completely screened pc
    //        p_upper = non-screened pc
    real_t Ecfree_term = Constants::ec * Ec_free[ir] /(Constants::me * Constants::c);
    real_t Ectot_term  = Constants::ec * Ec_tot[ir]  /(Constants::me * Constants::c);
    

    *p_lower = 1/sqrt(sqrt(pStar_params.constTerm)/Ecfree_term);

    // If pStar is smaller than p_lower (for some reason), reduce by 30%
    bool isPLoUnderestimate = (pStarFunction(*p_lower, &pStar_params) < 0);
    while(!isPLoUnderestimate){
        *p_upper = *p_lower;
        *p_lower *= 0.7;
        isPLoUnderestimate = (pStarFunction(*p_lower, &pStar_params) < 0);
    }
    if(!(*p_upper==0))
        return;

    *p_upper = 1/sqrt(sqrt(pStar_params.constTerm)/Ectot_term); 
    bool isPUpOverestimate = (pStarFunction(*p_upper, &pStar_params) > 0);
    while (!isPUpOverestimate){
        *p_upper *= 1.4;
        isPUpOverestimate = (pStarFunction(*p_upper, &pStar_params) > 0);
    }
    


    
}

// Takes a p interval [x_lower, x_upper] and iterates at most max_iter=30 times (or to a relative error of rel_error=0.001)
// to find an estimate for p_Star 
void CollisionQuantityHandler::FindPStarRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func){
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    gsl_root_fsolver_set (s, &gsl_func, x_lower, x_upper); 

    int status;
    real_t rel_error = 1e-3;
    len_t max_iter = 30;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status   = gsl_root_fsolver_iterate (s);
        *root    = gsl_root_fsolver_root (s);
        x_lower = gsl_root_fsolver_x_lower (s);
        x_upper = gsl_root_fsolver_x_upper (s);
        status   = gsl_root_test_interval (x_lower, x_upper, 0, rel_error);

        if (status == GSL_SUCCESS){
            gsl_root_fsolver_free(s);
            break;
        }
    }
}

/**
 * Calculates the runaway rate due to beta decay of tritium. (to be multiplied by n_tritium)
 */
real_t CollisionQuantityHandler::evaluateTritiumRate(len_t ir){
    if (criticalREMomentum[ir] == DBL_MAX)
        return 0;

    real_t tau_halfLife = 12.32 * 365.24 *24*60*60; // 12.32 years, in seconds

    real_t gamma_c = sqrt(1+criticalREMomentum[ir]*criticalREMomentum[ir]);

    real_t decayMaxEnergyEV = 18.6e3; // maximum beta electron kinetic energy 
    real_t w = Constants::mc2inEV * (gamma_c-1) / decayMaxEnergyEV;
    real_t fracAbovePc = 1 + sqrt(w)*( -(35/8)*w + (21/4)*w*w - (15/8)*w*w*w);

    return log(2) /* * n_tritium */ /tau_halfLife * fracAbovePc;
}


// Evaluates total cross section for Compton scattering into p>pc due to incident photon of energy Eg (units of mc and mc2)
// Eq (29) in Martin-Solis NF 2017
real_t CollisionQuantityHandler::evaluateComptonTotalCrossSectionAtP(real_t Eg, real_t pc){
    real_t x = Eg;
    real_t Wc = sqrt(1+pc*pc)-1;
    real_t cc = 1 - 1/Eg * Wc /( Eg - Wc );
    return M_PI * Constants::r0 * Constants::r0 * ( (x*x-2*x-2)/(x*x*x) * log( (1+2*x)/( 1+x*(1-cc) ) ) 
        + 1/(2*x) * ( 1/( (1+x*(1-cc))*(1+x*(1-cc)) ) - 1/( (1+2*x)*(1+2*x) ) ) 
        - 1/(x*x*x) * ( 1 - x - (1+2*x) / (1+x*(1-cc)) - x*cc )   );
}

// Photon spectral flux density, Eq (24) in Martin-Solis NF 2017
real_t CollisionQuantityHandler::evaluateComptonPhotonFluxSpectrum(real_t Eg){
    real_t ITERPhotonFluxDensity = 1e18; // 1/m^2s
    real_t z = (1.2 + log(Eg * Constants::mc2inEV/1e6) ) / 0.8;
    return ITERPhotonFluxDensity * exp( - exp(z) - z + 1 );
}


// The integrand in the evaluation of the total production rate integral(flux density * cross section ) 
struct ComptonParam {real_t pc; CollisionQuantityHandler *collQtyHand;};
real_t ComptonIntegrandFunc(real_t Eg, void *par){
    struct ComptonParam *params = (struct ComptonParam *) par;
    
    real_t pc = params->pc;
    CollisionQuantityHandler *collQtyHand = params->collQtyHand;

    return collQtyHand->evaluateComptonPhotonFluxSpectrum(Eg) * collQtyHand->evaluateComptonTotalCrossSectionAtP(Eg,pc);
}

// returns (dnRE/dt)_compton at radial index ir
real_t CollisionQuantityHandler::evaluateComptonRate(len_t ir,gsl_integration_workspace *gsl_ad_w){
    if(criticalREMomentum[ir]==DBL_MAX)
        return 0;

    struct ComptonParam  params= {criticalREMomentum[ir], this};
    gsl_function ComptonFunc;
    ComptonFunc.function = &(ComptonIntegrandFunc);
    ComptonFunc.params = &params;

    real_t gamma_c = sqrt(1+criticalREMomentum[ir]*criticalREMomentum[ir]);
    real_t Eg_min = (criticalREMomentum[ir] + gamma_c - 1) /2;
    real_t valIntegral;
    // qagiu assumes an infinite upper boundary
    real_t epsrel = 1e-4;
    real_t epsabs;
    gsl_integration_qagiu(&ComptonFunc, Eg_min , 0, epsrel, 1000, gsl_ad_w, &valIntegral, &epsabs);
    return n_tot[ir]*valIntegral;
}

// Evaluates and stores the avalanche, tritium and compton growth rates
void CollisionQuantityHandler::CalculateGrowthRates(){
    //len_t id_nRE = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
    //real_t *nRE = unknowns->GetUnknownData(id_nRE);


    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);


    len_t *tritiumIndices = nullptr;
    len_t numTritiumIndices = 0;

    real_t gamma_crit;
    avalancheRate = new real_t[this->nr];
    tritiumRate   = new real_t[this->nr];
    comptonRate   = new real_t[this->nr];
    for (len_t ir = 0; ir<this->nr; ir++){
        // we still haven't implemented the relativistic corrections in criticalREmomentum, 
        // but let's keep it like this for now in case we do in the future.

        const real_t nTritium = ionHandler->GetTritiumDensity(ir,tritiumIndices,numTritiumIndices);

        if(criticalREMomentum[ir]==DBL_MAX){
            avalancheRate[ir] = 0;
        } else {
            gamma_crit = sqrt( 1 + criticalREMomentum[ir]*criticalREMomentum[ir] );
            avalancheRate[ir] = 0.5 * n_tot[ir] * constPreFactor / (gamma_crit-1) ;
            //avalancheRate[ir] = 0.5 * nRE[ir] * constPreFactor / (gamma_crit-1) ;
        }
        tritiumRate[ir] = nTritium*evaluateTritiumRate(ir);
        comptonRate[ir] = evaluateComptonRate(ir,gsl_ad_w);
    }
}


real_t CollisionQuantityHandler::evaluateBremsStoppingForceAtP(len_t /*ir*/, real_t /*p*/){

    if(settings->bremsstrahlung_mode==OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER){
        throw NotImplementedException("Bremsstrahlung setting STOPPING_POWER not yet supported");
    } else if(settings->bremsstrahlung_mode==OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_BOLTZMANN){
        throw NotImplementedException("Bremsstrahlung setting BOLTZMANN not yet supported");
    } else {
        throw NotImplementedException("Bremsstrahlung setting not yet supported");
    }
    return 0;
}
