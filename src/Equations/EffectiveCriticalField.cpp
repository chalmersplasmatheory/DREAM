/**
 * Calculates the effective critical electric field
 * Settings parameters: 
 * - collQtySettings->collfreq_type (completely screened or not)
 * - Eceff_mode: 
 *  COLLQTY_ECEFF_MODE_EC_TOT      // Gives Ectot including all bound electrons (or Ec_free if no impurities/complete screening)
 *  COLLQTY_ECEFF_MODE_CYLINDRICAL // Sets Eceff using the Hesslow formula ignoring trapping effects.
 *  COLLQTY_ECEFF_MODE_SIMPLE      // An approximate numerical calculation with a simplified account of trapping effects
 *  COLLQTY_ECEFF_MODE_FULL        // Full 'Lehtinen theory' expression.
 * 
 * 
 */
#include "DREAM/Equations/EffectiveCriticalField.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
EffectiveCriticalField::EffectiveCriticalField(ParametersForEceff *par, AnalyticDistributionRE *analyticRE)
    : Eceff_mode(par->Eceff_mode), collSettingsForEc(par->collSettingsForEc), 
    rGrid(par->rGrid), nuS(par->nuS), nuD(par->nuD), ions(par->ions), lnLambda(par->lnLambda), 
    thresholdToNeglectTrappedContribution(par->thresholdToNeglectTrappedContribution), 
    fdfsolve(par->fdfsolve) 
{
    gsl_parameters.rGrid = par->rGrid;
    gsl_parameters.nuS = par->nuS;
    gsl_parameters.nuD = par->nuD;
    gsl_parameters.fgType = par->fgType;
    gsl_parameters.gsl_ad_w = par->gsl_ad_w;
    gsl_parameters.fmin = par->fmin;
    gsl_parameters.collSettingsForEc = par->collSettingsForEc;
//    gsl_parameters.QAG_KEY = GSL_INTEG_GAUSS31;
    gsl_parameters.analyticDist = analyticRE;
    gsl_parameters.BAFunc_par = nullptr;
}

/**
 * Destructor
 */
EffectiveCriticalField::~EffectiveCriticalField(){ 
    DeallocateQuantities();
}

/**
 * Deallocator
 */
void EffectiveCriticalField::DeallocateQuantities(){
    if(ECRIT_ECEFFOVERECTOT_PREV != nullptr){
        delete [] ECRIT_ECEFFOVERECTOT_PREV;
        delete [] ECRIT_POPTIMUM_PREV;

        if ((Eceff_mode == OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE) || (Eceff_mode == OptionConstants::COLLQTY_ECEFF_MODE_FULL)){
            for (len_t ir = 0; ir<nr; ir++) {
                delete [] EOverUnityContrib[ir]; 
                delete [] SynchOverUnityContrib[ir];

                gsl_spline_free (gsl_parameters.EContribSpline[ir]); 
                gsl_spline_free (gsl_parameters.SynchContribSpline[ir]); 
                gsl_interp_accel_free (gsl_parameters.EContribAcc[ir]);
                gsl_interp_accel_free (gsl_parameters.SynchContribAcc[ir]);
            }

            delete [] EOverUnityContrib;
            delete [] SynchOverUnityContrib;

            delete [] gsl_parameters.EContribSpline;
            delete [] gsl_parameters.SynchContribSpline;

        }
    }
}

/** 
 * To be called when the grid has been rebuilt; 
 * will reallocate memory for and calculates 
 * grid-dependent quantities
 */
bool EffectiveCriticalField::GridRebuilt(){
    DeallocateQuantities();
    nr = rGrid->GetNr(); // update afterwards so gsl_free doesn't crash
    ECRIT_ECEFFOVERECTOT_PREV = new real_t[nr];
    ECRIT_POPTIMUM_PREV = new real_t[nr];
    // Initial guess: Eceff/Ectot \approx 1.0, 
    // with a corresponding critical momentum at p=10mc    
    for(len_t ir=0; ir<nr; ir++){ 
        ECRIT_ECEFFOVERECTOT_PREV[ir] = 1.0;
        ECRIT_POPTIMUM_PREV[ir] = 10;
    }

    gsl_parameters.Eterm = 0;
    gsl_parameters.p = 0;
    gsl_parameters.p_ex_lo = 0;
    gsl_parameters.p_ex_up = 0; 


    if ((Eceff_mode == OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE) || (Eceff_mode == OptionConstants::COLLQTY_ECEFF_MODE_FULL)){
        this->EOverUnityContrib = new real_t*[nr];
        this->SynchOverUnityContrib = new real_t*[nr];

        this->gsl_parameters.EContribSpline = new gsl_spline*[nr]; 
        this->gsl_parameters.SynchContribSpline = new gsl_spline*[nr];

        this->gsl_parameters.EContribAcc =  new gsl_interp_accel*[nr]; // the accelerators cache values from the splines
        this->gsl_parameters.SynchContribAcc = new gsl_interp_accel*[nr]; 

        real_t fracPointsLower   = 0.4; 
        real_t fracUpperInterval = 1.0/3.0;

        /**
         * Create X_vec array on [0,1], on which we will interpolate 
         * the RE pitch distribution averaged coefficients.
         * Constructed as:
         *  N_A_VALUES * fracPointsLower on the interval [0,1-fracUpperInterval]
         *  remaining points on the interval [1-fracUpperInterval, 1]
         * corresponding (here) to a denser grid near X_vec ~ 1, corresponding to
         * large A (strong electric fields, beam-like distributions)
         */
        len_t N1 = (len_t) (N_A_VALUES * fracPointsLower); // rounds down to a natural number
        len_t N2 = N_A_VALUES - N1;
        for(len_t i=0; i<N1; i++)
            X_vec[i] = i*(1-fracUpperInterval)/(N1-1);
        for(len_t i=1; i<=N2; i++)
            X_vec[N1-1+i] = 1 - fracUpperInterval + i*fracUpperInterval/N2;

        real_t synchrotronPrefactor = Constants::ec * Constants::ec * Constants::ec * Constants::ec 
                                / ( 6 * M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me
                                * Constants::c * Constants::c * Constants::c);
        for (len_t ir = 0; ir<nr; ir++){
            EOverUnityContrib[ir] = new real_t[N_A_VALUES];
            SynchOverUnityContrib[ir] = new real_t[N_A_VALUES];
            gsl_parameters.ir = ir;
            real_t FSA_B, FSA_B2, Bmin; // get grid parameters evaluated on the appropriate radial grid
            if(gsl_parameters.fgType==FVM::FLUXGRIDTYPE_RADIAL){
                FSA_B = rGrid->GetFSA_B_f(ir);
                FSA_B2 = rGrid->GetFSA_B2_f(ir);
                Bmin = rGrid->GetBmin_f(ir);
            } else {
                FSA_B = rGrid->GetFSA_B(ir);
                FSA_B2 = rGrid->GetFSA_B2(ir);
                Bmin = rGrid->GetBmin(ir);
            }
            // store constants that do not depend on plasma parameters or p
            gsl_parameters.CONST_E = Constants::ec  / (Constants::me * Constants::c) * sqrt(FSA_B2);
            gsl_parameters.CONST_EFact = gsl_parameters.CONST_E / FSA_B;
            gsl_parameters.CONST_Synch = synchrotronPrefactor * Bmin*Bmin;
            for (len_t i = 0; i<N_A_VALUES-1; i++){
                gsl_parameters.A = GetAFromX(X_vec[i]);               
                CreateLookUpTableForUIntegrals(&gsl_parameters, EOverUnityContrib[ir][i], SynchOverUnityContrib[ir][i]);
            }
            // known values at A=inf (X=1)
            EOverUnityContrib[ir][N_A_VALUES-1] = 1;
            SynchOverUnityContrib[ir][N_A_VALUES-1] = 0;
            
            gsl_parameters.EContribAcc[ir] = gsl_interp_accel_alloc();
            gsl_parameters.EContribSpline[ir] = gsl_spline_alloc (gsl_interp_steffen, N_A_VALUES);
            gsl_spline_init (gsl_parameters.EContribSpline[ir], X_vec, EOverUnityContrib[ir], N_A_VALUES);

            gsl_parameters.SynchContribAcc[ir] = gsl_interp_accel_alloc();
            gsl_parameters.SynchContribSpline[ir] = gsl_spline_alloc (gsl_interp_steffen, N_A_VALUES);
            gsl_spline_init (gsl_parameters.SynchContribSpline[ir], X_vec, SynchOverUnityContrib[ir], N_A_VALUES);

        }
    }
    return true;
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
    switch (Eceff_mode)
    {
        case OptionConstants::COLLQTY_ECEFF_MODE_EC_TOT : { // or COLLQTY_ECEFF_MODE_NOSCREENING to be consistent with 
                                                            // for example GetConnorHastieField_NOSCREENING?
            if(collSettingsForEc->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED)
                for(len_t ir=0; ir<nr; ir++)
                    effectiveCriticalField[ir] = Ec_free[ir];
            else
                for(len_t ir=0; ir<nr; ir++)
                    effectiveCriticalField[ir] = Ec_tot[ir]; 
        } 
        break;
        case OptionConstants::COLLQTY_ECEFF_MODE_CYLINDRICAL : 
            for(len_t ir=0; ir<nr; ir++)
                effectiveCriticalField[ir] = CalculateEceffPPCFPaper(ir);
        break;
        case OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE : 
            [[fallthrough]];
        case OptionConstants::COLLQTY_ECEFF_MODE_FULL : {  
            gsl_function_fdf UExtremumFunc;
            for (len_t ir=0; ir<nr; ir++){
                gsl_parameters.ir = ir;
                // it was found empirically that with a 10% margin, typical simulations
                // will seldom end up outside of the interval
                gsl_parameters.p_ex_lo = 0.9 * ECRIT_POPTIMUM_PREV[ir];
                gsl_parameters.p_ex_up = 1.1 * ECRIT_POPTIMUM_PREV[ir];

                UExtremumFunc.f = &(FindUExtremumAtE);
                UExtremumFunc.df = &(FindUExtremumAtE_df);
                UExtremumFunc.fdf = &(FindUExtremumAtE_fdf);
                UExtremumFunc.params = &gsl_parameters; 

                real_t E_root = ECRIT_ECEFFOVERECTOT_PREV[ir] * Ec_tot[ir];
                RunawayFluid::FindRoot_fdf(E_root, UExtremumFunc,fdfsolve);
                effectiveCriticalField[ir] = E_root;
                ECRIT_ECEFFOVERECTOT_PREV[ir] = effectiveCriticalField[ir]/Ec_tot[ir];
                ECRIT_POPTIMUM_PREV[ir] = gsl_parameters.p_optimum;                
            }
        }
        break;
        default :
            break;
    }
}

/**
 * Calculates the effective critical field using Eqs (23)-(24) in Hesslow et al, PPCF 60, 074010 (2018),
 * which is also the formula published on GitHub.
 */
real_t EffectiveCriticalField::CalculateEceffPPCFPaper(len_t ir){
    real_t  lnLambdaC = lnLambda->GetLnLambdaC(ir);
    real_t  ne_free = ions->GetFreeElectronDensityFromQuasiNeutrality(ir); 
    real_t  ne_tot = ions->GetFreePlusBoundElectronDensity(ir);
    real_t  Zeff = ions->GetZeff(ir); 
    real_t  Zfulleff = 0; 
    real_t  Ne2_nj = 0;
    real_t  Ne_nj = 0; 

    // Depending on settings, Ectot can be with the thermal value so we re-calculate here just in case.
    // I think it should always be lnLc here. 
    real_t Ec_free_lnLambdaC = ne_free * lnLambdaC * 
        4*M_PI *Constants::r0*Constants::r0 *Constants::me * Constants::c*Constants::c / Constants::ec; 

    len_t ionOffset = 0;
    len_t Z;
    real_t nj_Over_ne_free;
    for (len_t iz = 0; iz < ions->GetNZ(); iz++){ 
        Z = ions->GetZ(iz);
        for (len_t Z0 = 0; Z0 <= Z; Z0++, ionOffset++){ 
            nj_Over_ne_free = ions->GetIonDensity(ir, iz, Z0) / ne_free;
            Zfulleff += Z*Z*nj_Over_ne_free;
            Ne2_nj += (Z-Z0)*(Z-Z0)*nj_Over_ne_free;
            Ne_nj += (Z-Z0)*nj_Over_ne_free;
        }
    }
    
    real_t nuD0 = (1 + Zeff) - 2.0/3.0 * Ne2_nj /lnLambdaC;
    real_t nuD1 = Zfulleff/lnLambdaC;
    real_t nuS0 = 1;
    real_t nuS1 = (1 + 3.0 * Ne_nj)/(2*lnLambdaC);
    ionOffset = 0;
    for (len_t iz = 0; iz < ions->GetNZ(); iz++){ 
        Z = ions->GetZ(iz);
        for (len_t Z0 = 0; Z0 < Z; Z0++, ionOffset++){ // exclude Z0=Z since they are not partially ionized
            real_t ajBar = nuD->GetIonEffectiveSizeAj(iz, Z0);
            real_t IjBar = nuS->GetMeanExcitationEnergy(iz, Z0);
            nj_Over_ne_free = ions->GetIonDensity(ir, iz, Z0) /  ne_free;
            nuD0 += (Z*Z - Z0*Z0) * log(ajBar) * nj_Over_ne_free / lnLambdaC;
            nuS0 += (Z-Z0) * (-log(IjBar)-1) * nj_Over_ne_free / lnLambdaC;
        }
    }
    // synchrotron radiation term
    real_t B2 = rGrid->GetFSA_B2(ir) * rGrid->GetBmin(ir) * rGrid->GetBmin(ir); 
    real_t tauRadInv =  B2/(ne_free/1e20)/(15.44*lnLambdaC);
    
    // bremsstrahlung term
    real_t bremsprefactor = Constants::alpha * Zfulleff/lnLambdaC;
    // partially ionized expression based on an approximate bremsstrahlung formula used in the paper, but we later found out that 
    // the introduced error in the approximation was larger than the screening effects.
    // Keep this here for now for comparison (agrees to at least the first five decimals with Matlab scripts)
    // real_t phib1 = bremsprefactor*0.35;
    // real_t phib2 = bremsprefactor*0.20;

    // non-screened limit. The stopping-power comes from Kock&Motz: \sum_i p*c*ni * [4BN(b)] 
    real_t phib1 = bremsprefactor*(log(2)-1/3)/M_PI;
    real_t phib2 = bremsprefactor/M_PI;

    real_t Eceff = ne_tot/ne_free; // Eceff0
    real_t delta = 0;
    real_t pc0 = nuD0/(2*nuS1);

    auto calcEceff = [&] () { Eceff =  nuS0 + nuS1 *( (1+nuD1/nuD0)*log(pc0) + sqrt(2*delta+1)); };
    auto calcDelta = [&] () { delta = (nuD0/(nuS1*nuS1)) * (nuD0*tauRadInv/Eceff+phib1+phib2*log(pc0)); };

    calcDelta(); // delta0
    calcEceff(); // Eceff1
    calcDelta(); // delta1
    calcEceff(); // Eceff
    return Ec_free_lnLambdaC * Eceff;
}

/**
 * Returns the minimum of the acceleration function -U (with respect to p) at a given Eterm 
 */
real_t EffectiveCriticalField::FindUExtremumAtE(real_t Eterm, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    params->Eterm = Eterm;
    gsl_min_fminimizer *gsl_fmin = params->fmin;
    real_t 
        p_ex_guess, p_ex_lo, p_ex_up,
        F_ex_guess, F_ex_lo, F_ex_up; // function values

    gsl_function F;
    F.function = &(UAtPFunc);
    F.params = params;
    real_t p_upper_threshold = 1000; // larger momenta are not physically relevant in our scenarios
    FindPExInterval(
        p_ex_guess, p_ex_lo, p_ex_up,
        F_ex_guess, F_ex_lo, F_ex_up,
        p_upper_threshold,params
    );

    // If the extremum is at a larger momentum than p_upper_threshold (or doesn't exist at all), 
    // we will define Eceff as the value where U(p_upper_threshold) = 0. 
    if(p_ex_up > p_upper_threshold)
        return UAtPFunc(p_upper_threshold,params);

    gsl_min_fminimizer_set_with_values(
        gsl_fmin, &F, 
        p_ex_guess, F_ex_guess, 
        p_ex_lo, F_ex_lo, 
        p_ex_up, F_ex_up
    );

    int status;
    real_t rel_error = 1e-2, abs_error=0;
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
    params->p_optimum = p_ex_guess;
    real_t minimumFValue = gsl_min_fminimizer_f_minimum(gsl_fmin);
    return minimumFValue;
}
/**
 * Returns a finite-difference differentiation of the 'FindUExtremumAtE' function
 */
real_t EffectiveCriticalField::FindUExtremumAtE_df(real_t Eterm, void *par){
    real_t h = 0.01*Eterm;
    return (FindUExtremumAtE(Eterm,par) - FindUExtremumAtE(Eterm-h,par))/h;
}
/**
 * Stores a finite-difference differentiation of the 'FindUExtremumAtE' function
 * as 'df' and the function value 'FindUExtremumAtE' as 'f'
 */
void EffectiveCriticalField::FindUExtremumAtE_fdf(real_t Eterm, void *par, real_t *f, real_t *df){
    real_t h = 0.01*Eterm;
    *f = FindUExtremumAtE(Eterm,par);
    *df = (*f-FindUExtremumAtE(Eterm-h,par)) / h;
}

/**
 * Finds an interval p \in [p_ex_lower, p_ex_upper] in which a minimum of -U(p) exists.  
 */
void EffectiveCriticalField::FindPExInterval(
    real_t &p_ex_guess, real_t &p_ex_lower, real_t &p_ex_upper, 
    real_t &F_ex_guess, real_t &F_ex_lower, real_t &F_ex_upper, 
    real_t p_upper_threshold, UContributionParams *params
){
    p_ex_lower = params->p_ex_lo;
    p_ex_upper = params->p_ex_up;
    p_ex_guess = sqrt(p_ex_lower * p_ex_upper);
    F_ex_lower = UAtPFunc(p_ex_lower,params);
    F_ex_upper = UAtPFunc(p_ex_upper,params);
    F_ex_guess = UAtPFunc(p_ex_guess,params);
    
    if( (F_ex_guess < F_ex_upper) && (F_ex_guess < F_ex_lower) ) // at least one minimum exists on the interval
        return;
    else if ( F_ex_guess > F_ex_lower) // Minimum located at p<p_ex_guess
        while(F_ex_guess > F_ex_lower){
            p_ex_upper = p_ex_guess;
            p_ex_guess = p_ex_lower;
            p_ex_lower /= 2;
            F_ex_guess = F_ex_lower; //UAtPFunc(*p_ex_guess,params);
            F_ex_lower = UAtPFunc(p_ex_lower,params);
        }
    else // Minimum at p>p_ex_guess
        while( (F_ex_guess > F_ex_upper) && (p_ex_upper < p_upper_threshold)){
            p_ex_lower = p_ex_guess;
            p_ex_guess = p_ex_upper;
            p_ex_upper *= 2;
            F_ex_guess = F_ex_upper;//UAtPFunc(*p_ex_guess,params);
            F_ex_upper = UAtPFunc(p_ex_upper,params);
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
real_t UPartialContributionForInterpolation(real_t xi0, void *par){
    struct EffectiveCriticalField::UContributionParams *params = (struct EffectiveCriticalField::UContributionParams *) par;
    len_t ir = params->ir;
    return params->preFactorFunc(xi0)*params->rGrid->EvaluatePXiBounceIntegralAtP(ir,xi0,params->fgType,params->BAFunc,params->BAFunc_par,params->BAList)
        * params->analyticDist->evaluatePitchDistributionFromA(ir,xi0,params->A);
}

void EffectiveCriticalField::CreateLookUpTableForUIntegrals(UContributionParams *params, real_t &EContrib, real_t &SynchContrib){
    FVM::RadialGrid *rGrid = params->rGrid;
    len_t ir = params->ir;
    FVM::fluxGridType fluxGridType = params->fgType;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    
    
    real_t xiT = (fluxGridType == FVM::FLUXGRIDTYPE_RADIAL) ?
            rGrid->GetXi0TrappedBoundary_fr(ir) :
            rGrid->GetXi0TrappedBoundary(ir);
    if(xiT < thresholdToNeglectTrappedContribution)
        xiT = 0;
    
    // Evaluates the contribution from electric field term A^p coefficient
    params->BAFunc = FVM::RadialGrid::BA_FUNC_XI;
    params->BAList = FVM::RadialGrid::BA_PARAM_XI;
    params->preFactorFunc = [](real_t xi0){return xi0;};
    real_t error;
    real_t epsabs = 1e-6, epsrel = 5e-3, lim = gsl_ad_w->limit; 
    gsl_function GSL_func;
    GSL_func.function = &(UPartialContributionForInterpolation);
    GSL_func.params = params;

    // size_t key = GSL_INTEG_GAUSS31; // integration order w qag in interval 1-6, where 6 is the highest order
    // have tried three options here: 
    // - qags (works, but probably slower than necessary)
    // - qagp (I get SEGABRT, don't understand why)
    // - qag, should work for xiT = 0 since there are no singularities then. chose key between 1 and 6
    //      GSL_INTEG_GAUSS11 => error, divide by zero
    //      GSL_INTEG_GAUSS21 => almost the same as qags
    //      GSL_INTEG_GAUSS31 => a little slower than qags
    //      GSL_INTEG_GAUSS51 => factor 2 slower than qags
    // note that the RunawayFluid currently doesn't test trapping effects, so the xiT>0 cases are not evaluated
    if(xiT){
        real_t EContrib1, EContrib2;
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&EContrib1,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&EContrib2,&error);
        EContrib = EContrib1 + EContrib2;
    } else
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&EContrib,&error); 
    // Evaluates the contribution from slowing down term A^p coefficient
    params->BAFunc = FVM::RadialGrid::BA_FUNC_UNITY;
    params->BAList = FVM::RadialGrid::BA_PARAM_UNITY;
    params->preFactorFunc = [](real_t){return 1;};
    real_t UnityContrib;    
    if(xiT){
        real_t UnityContrib1, UnityContrib2, UnityContrib3;
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib1,&error);
        gsl_integration_qags(&GSL_func,0,xiT,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib2,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib3,&error);
        UnityContrib = UnityContrib1 + UnityContrib2 + UnityContrib3;
    } else
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib,&error);

    // Evaluates the contribution from synchrotron term A^p coefficient
    params->BAFunc = FVM::RadialGrid::BA_FUNC_B_CUBED;
    params->BAList = FVM::RadialGrid::BA_PARAM_B_CUBED;
    params->preFactorFunc = [](real_t xi0){return 1-xi0*xi0;};
    if(xiT){
        real_t SynchContrib1, SynchContrib2, SynchContrib3;
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib1,&error);
        gsl_integration_qags(&GSL_func,0,xiT,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib2,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib3,&error);
        SynchContrib = SynchContrib1 + SynchContrib2 + SynchContrib3;
    } else
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib,&error);
    EContrib /= UnityContrib;
    SynchContrib /= UnityContrib;
}


/**
 * Evaluates -U(p) at given Eterm.
 */
real_t EffectiveCriticalField::UAtPFunc(real_t p, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    params->p = p;
    len_t ir = params->ir;
    real_t Eterm = params->Eterm;
    SlowingDownFrequency *nuS = params->nuS;
    PitchScatterFrequency *nuD = params->nuD;
    CollisionQuantity::collqty_settings *collSettingsForEc = params->collSettingsForEc;    

    real_t E = params->CONST_E * Eterm; 
    real_t pNuD = p*nuD->evaluateAtP(ir,p,collSettingsForEc);    
    real_t A = 2*E/pNuD;
    
    // Evaluates the contribution from electric field term A^p coefficient
    real_t Efactor = params->CONST_EFact * Eterm; 
    real_t SynchrotronFactor = -p*sqrt(1+p*p) * params->CONST_Synch; 

    real_t EContrib = Efactor * gsl_spline_eval(params->EContribSpline[ir], GetXFromA(A), params->EContribAcc[ir]);
    real_t NuSContrib = -p*nuS->evaluateAtP(ir,p,collSettingsForEc);
    real_t SynchContrib = SynchrotronFactor * gsl_spline_eval(params->SynchContribSpline[ir], GetXFromA(A), params->SynchContribAcc[ir]);

    return -(EContrib + NuSContrib + SynchContrib) ;
}
