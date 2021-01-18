/**
 * Implementation of a class that calculates and stores quantities related to the runaway growth rate, 
 * e.g. the effective critical field and runaway momentum, as well as avalanche and Dreicer growths etc.
 */

#include <string>
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
#include "DREAM/Equations/ConnorHastie.hpp"
#include "DREAM/Equations/DreicerNeuralNetwork.hpp"
#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/TimeKeeper.hpp"
#include <limits>

using namespace DREAM;

const real_t RunawayFluid::tritiumHalfLife =  3.888e8;    // 12.32 years, in seconds
const real_t RunawayFluid::tritiumDecayEnergyEV = 18.6e3; // maximum beta electron kinetic energy in eV 

const len_t  RunawayFluid::conductivityLenT = 14;
const len_t  RunawayFluid::conductivityLenZ = 6;
//const real_t RunawayFluid::conductivityBraams[conductivityLenZ*conductivityLenT] = {3.75994, 3.7549, 3.7492, 3.72852, 3.6842, 3.57129, 3.18206, 2.65006, 2.03127, 1.33009, 0.94648, 0.67042, 0.42422, 0.29999, 7.42898, 7.27359, 7.12772, 6.73805, 6.20946, 5.43667, 4.13733, 3.13472, 2.27862, 1.45375, 1.02875, 0.72743, 0.46003, 0.32528, 8.7546, 8.53281, 8.32655, 7.78445, 7.06892, 6.06243, 4.47244, 3.32611, 2.39205, 1.51805, 1.07308, 0.75853, 0.47965, 0.33915, 10.39122, 10.07781, 9.78962, 9.04621, 8.09361, 6.80431, 4.8805, 3.57303, 2.54842, 1.61157, 1.13856, 0.80472, 0.50885, 0.35979, 11.33006, 10.95869, 10.61952, 9.75405, 8.66306, 7.21564, 5.11377, 3.72206, 2.64827, 1.67382, 1.18263, 0.83593, 0.52861, 0.37377, 12.76615, 12.29716, 11.87371, 10.81201, 9.50746, 7.82693, 5.47602, 3.96944, 2.82473, 1.7887, 1.2649, 0.89443, 0.56569, 0.4};
const real_t RunawayFluid::conductivityBraams[conductivityLenZ*conductivityLenT] = {12.7661,12.2972,11.8737,10.8120,9.5075,7.8269,5.4760,3.9694,2.8247,1.7887,1.2649,0.8944,0.5657,0.4000,11.3301,10.9587,10.6195,9.7540,8.6631,7.2156,5.1138,3.7221,2.6483,1.6738,1.1826,0.8359,0.5286,0.3738,10.3912,10.0778,9.7896,9.0462,8.0936,6.8043,4.8805,3.5730,2.5484,1.6116,1.1386,0.8047,0.5089,0.3598,8.7546,8.5328,8.3265,7.7844,7.0689,6.0624,4.4724,3.3261,2.3920,1.5180,1.0731,0.7585,0.4797,0.3392,7.4290,7.2736,7.1277,6.7381,6.2095,5.4367,4.1373,3.1347,2.2786,1.4538,1.0288,0.7274,0.4600,0.3253,3.7599,3.7549,3.7492,3.7285,3.6842,3.5713,3.1821,2.6501,2.0313,1.3301,0.9465,0.6704,0.4242,0.3000};
const real_t RunawayFluid::conductivityTmc2[conductivityLenT] = {0,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100};
const real_t RunawayFluid::conductivityX[conductivityLenZ]    = {0,0.090909090909091,0.166666666666667,0.333333333333333,0.5,1};

/**
 * Constructor.
 */
RunawayFluid::RunawayFluid(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, SlowingDownFrequency *nuS, 
    PitchScatterFrequency *nuD, CoulombLogarithm *lnLee,
    CoulombLogarithm *lnLei, CollisionQuantity::collqty_settings *cqs,
    IonHandler *ions, AnalyticDistributionRE *distRE,
    OptionConstants::conductivity_mode cond_mode,
    OptionConstants::eqterm_dreicer_mode dreicer_mode,
    OptionConstants::collqty_Eceff_mode Eceff_mode,
    OptionConstants::eqterm_avalanche_mode ava_mode,
    OptionConstants::eqterm_compton_mode compton_mode,
    real_t compton_photon_flux
) : nuS(nuS), nuD(nuD), lnLambdaEE(lnLee), lnLambdaEI(lnLei),
    unknowns(u), ions(ions), analyticRE(distRE), cond_mode(cond_mode), dreicer_mode(dreicer_mode),
    Eceff_mode(Eceff_mode), ava_mode(ava_mode), compton_mode(compton_mode),
    compton_photon_flux(compton_photon_flux)
 {
    this->gridRebuilt = true;
    this->rGrid = g->GetRadialGrid();

    id_ncold = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ntot  = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_ni    = this->unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Tcold = this->unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Eterm = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    id_jtot  = this->unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);

    this->gsl_ad_w = gsl_integration_workspace_alloc(1000);
    this->gsl_ad_w2 = gsl_integration_workspace_alloc(1000);
    this->fsolve = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    this->fdfsolve = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant);

    this->fmin = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

    // Set collision settings for the Eceff calculation; always include bremsstrahlung and
    // energy-dependent Coulomb logarithm. The user only chooses collfreq_type, in practice.
    collSettingsForEc = new CollisionQuantity::collqty_settings;
    collSettingsForEc->collfreq_mode       = cqs->collfreq_mode;
    collSettingsForEc->collfreq_type       = cqs->collfreq_type;
    collSettingsForEc->pstar_mode          = cqs->pstar_mode;
    collSettingsForEc->screened_diffusion  = cqs->screened_diffusion;
    collSettingsForEc->lnL_type            = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    collSettingsForEc->bremsstrahlung_mode = OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;

    real_t thresholdToNeglectTrapped = 100*sqrt(std::numeric_limits<real_t>::epsilon());
    EffectiveCriticalField::ParametersForEceff par = {
        rGrid, nuS, nuD, FVM::FLUXGRIDTYPE_DISTRIBUTION, gsl_ad_w, gsl_ad_w2, fmin, collSettingsForEc,
        fdfsolve, Eceff_mode,ions,lnLambdaEI,thresholdToNeglectTrapped
    };
    this->effectiveCriticalFieldObject = new EffectiveCriticalField(&par, analyticRE);

    // Set collision settings for the critical-momentum calculation: takes input settings but 
    // ignores bremsstrahlung since generation never occurs at ultrarelativistic energies where it matters
    collSettingsForPc = new CollisionQuantity::collqty_settings;
    collSettingsForPc->collfreq_mode       = cqs->collfreq_mode;
    collSettingsForPc->collfreq_type       = cqs->collfreq_type;
    collSettingsForPc->pstar_mode          = cqs->pstar_mode;
    collSettingsForPc->screened_diffusion  = cqs->screened_diffusion;
    collSettingsForPc->lnL_type            = cqs->lnL_type;
    collSettingsForPc->bremsstrahlung_mode = OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT;

    // We always construct a Connor-Hastie runaway rate object, even if
    // the user would prefer to use the neural network. This is so that
    // we have something to fall back to in case the simulation enters
    // into a domain in which the neural network is not applicable.
    // (The ConnorHastie object is very cheap; it occupies a little more
    // than '2*nr' doubles and incurs no extra computational overhead)
    this->dreicer_ConnorHastie = new ConnorHastie(this);
    this->dreicer_ConnorHastie->IncludeCorrections(dreicer_mode == OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE);

    // Only construct neural network if explicitly requested...
    if (dreicer_mode == OptionConstants::EQTERM_DREICER_MODE_NEURAL_NETWORK
     || dreicer_mode == OptionConstants::EQTERM_DREICER_MODE_NONE)
        this->dreicer_nn = new DreicerNeuralNetwork(this);

    const gsl_interp2d_type *gsl_T = gsl_interp2d_bilinear; 
    gsl_cond = gsl_interp2d_alloc(gsl_T, conductivityLenT,conductivityLenZ);
    gsl_xacc = gsl_interp_accel_alloc();
    gsl_yacc = gsl_interp_accel_alloc();

    gsl_interp2d_init(gsl_cond, conductivityTmc2, conductivityX, conductivityBraams,conductivityLenT,conductivityLenZ);

    // Initialize timers
    this->timeKeeper = new FVM::TimeKeeper("RunawayFluid");
    this->timerTot = this->timeKeeper->AddTimer("total", "Total time");
    this->timerLnLambdaEE = this->timeKeeper->AddTimer("lnlambdaEE", "e-e Coulomb logarithm");
    this->timerLnLambdaEI = this->timeKeeper->AddTimer("lnlambdaEI", "e-i Coulomb logarithm");
    this->timerNuS = this->timeKeeper->AddTimer("nus", "Slowing-down frequency");
    this->timerNuD = this->timeKeeper->AddTimer("nud", "Pitch scattering frequency");
    this->timerDerived = this->timeKeeper->AddTimer("derived", "Derived quantities");
    this->timerEcEff = this->timeKeeper->AddTimer("eceff", "Effective critical electric field");
    this->timerPCrit = this->timeKeeper->AddTimer("pcrit", "Critical momentum");
    this->timerGrowthrates = this->timeKeeper->AddTimer("growthrates", "Runaway growthrates");    
}


/**
 * Destructor.
 */
RunawayFluid::~RunawayFluid(){
    DeallocateQuantities();

    gsl_integration_workspace_free(gsl_ad_w);
    gsl_integration_workspace_free(gsl_ad_w2);
    gsl_root_fsolver_free(fsolve);
    gsl_root_fdfsolver_free(fdfsolve);
    gsl_min_fminimizer_free(fmin);

    gsl_interp2d_free(gsl_cond);
    gsl_interp_accel_free(gsl_xacc);
    gsl_interp_accel_free(gsl_yacc);
    
    if (dreicer_ConnorHastie != nullptr)
        delete dreicer_ConnorHastie;
    if (dreicer_nn != nullptr)
        delete dreicer_nn;

    delete effectiveCriticalFieldObject;

    delete collSettingsForEc;
    delete collSettingsForPc;

    delete timeKeeper;
}


/**
 * Rebuilds all runaway quantities if plasma parameters have changed.
 */
void RunawayFluid::Rebuild(){
    this->timeKeeper->StartTimer(timerTot);

    // Macro for running accumulating timers
    #define TIME(NAME, STM) \
        do { this->timeKeeper->StartTimer( timer ## NAME ); (STM); this->timeKeeper->StopTimer( timer ## NAME ); } while (false)

    if(!parametersHaveChanged())
        return;
        
    if(gridRebuilt){
        nr = rGrid->GetNr();
        AllocateQuantities();
        effectiveCriticalFieldObject->GridRebuilt();
        gridRebuilt = false;
    }
    ncold   = unknowns->GetUnknownData(id_ncold);
    ntot    = unknowns->GetUnknownData(id_ntot);
    Tcold   = unknowns->GetUnknownData(id_Tcold);
    Eterm   = unknowns->GetUnknownData(id_Eterm);
    
    // The collision frequencies (nuS, nuD) uses the coulomb logarithms 
    // in a way that requires them to be rebuilt here.
    TIME(LnLambdaEE, lnLambdaEE->RebuildRadialTerms());
    TIME(LnLambdaEI, lnLambdaEI->RebuildRadialTerms());
    TIME(NuS, nuS->RebuildRadialTerms());
    TIME(NuD, nuD->RebuildRadialTerms());

    TIME(Derived, CalculateDerivedQuantities());
    TIME(EcEff, effectiveCriticalFieldObject->CalculateEffectiveCriticalField(Ec_tot, Ec_free,effectiveCriticalField));
    TIME(PCrit, CalculateCriticalMomentum());
    TIME(Growthrates, CalculateGrowthRates());

    this->timeKeeper->StopTimer(timerTot);
}


/** 
 * Returns true if any unknown quantities that affect runaway rates have changed. 
 */
bool RunawayFluid::parametersHaveChanged(){
    return unknowns->HasChanged(id_ncold) || unknowns->HasChanged(id_Tcold) || unknowns->HasChanged(id_ni) || 
            unknowns->HasChanged(id_Eterm) || gridRebuilt;
}

/**
 * Calculates the Connor-Hastie field Ec using the relativistic lnLambda
 * and the Dreicer field ED using the thermal lnLambda.
 */
void RunawayFluid::CalculateDerivedQuantities(){
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    for (len_t ir=0; ir<nr; ir++){
        real_t lnLc = lnLambdaEE->evaluateLnLambdaC(ir);
        real_t lnLT = lnLambdaEE->evaluateLnLambdaT(ir);
        // if running with lnLambda = THERMAL, override the relativistic lnLambda
        if(collSettingsForPc->lnL_type == OptionConstants::COLLQTY_LNLAMBDA_THERMAL)
            lnLc = lnLT;
        Ec_free[ir] = lnLc * ncold[ir] * constPreFactor * Constants::me * Constants::c / Constants::ec;
        Ec_tot[ir]  = lnLc * ntot[ir]  * constPreFactor * Constants::me * Constants::c / Constants::ec;
        EDreic[ir]  = lnLT * ncold[ir] * constPreFactor * Constants::me * Constants::c / Constants::ec * (Constants::mc2inEV / T_cold[ir]);
        
        if(ncold[ir] > 0){
            tauEERel[ir] = 1/(lnLc * ncold[ir] * constPreFactor); // = m*c/(e*Ec_free)
            tauEETh[ir]  = 1/(lnLT * ncold[ir] * constPreFactor) * pow(2*T_cold[ir]/Constants::mc2inEV,1.5); 
        } else { // if ncold=0 (for example at t=0 of hot tail simulation), set to infinite
            tauEERel[ir] = std::numeric_limits<real_t>::infinity();
            tauEETh[ir]  = std::numeric_limits<real_t>::infinity();
        }

        switch(cond_mode) {
            case OptionConstants::CONDUCTIVITY_MODE_BRAAMS: 
                electricConductivity[ir] = evaluateBraamsElectricConductivity(ir);
                break;
            case OptionConstants::CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS: 
                electricConductivity[ir] = evaluateSauterElectricConductivity(ir, true);
                break;
            case OptionConstants::CONDUCTIVITY_MODE_SAUTER_COLLISIONAL: 
                electricConductivity[ir] = evaluateSauterElectricConductivity(ir, false);
                break;
            default:
                break; // throw exception?            
        }
    }
}

/**
 * Is called to specify that the grid has been rebuilt, and that reallocation of all stored quantities is needed.
 */
void RunawayFluid::GridRebuilt(){
    gridRebuilt = true;
    lnLambdaEE->GridRebuilt();
    lnLambdaEI->GridRebuilt();
    nuS->GridRebuilt();
    nuD->GridRebuilt();    
    effectiveCriticalFieldObject->GridRebuilt();
}

/**
 * Finds the root of the provided gsl_function in the interval x_lower < root < x_upper. 
 * Is used both in the Eceff and pCrit calculations. 
 */
void RunawayFluid::FindRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func, gsl_root_fsolver *s){
    gsl_root_fsolver_set(s, &gsl_func, x_lower, x_upper); 
    int status;
    real_t epsrel = 1e-3;
    len_t max_iter = 30;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        gsl_root_fsolver_iterate (s);
        *root   = gsl_root_fsolver_root (s);
        x_lower = gsl_root_fsolver_x_lower (s);
        x_upper = gsl_root_fsolver_x_upper (s);
        status  = gsl_root_test_interval (x_lower, x_upper, 0, epsrel);
        if (status == GSL_SUCCESS)
            break;
    }
}

/**
 * Finds the root of the provided gsl_function_fdf using the provided
 * derivative-based solver.
 *  root: guess for the solution (and is overwritten by the obtained numerical solution)
 */
void RunawayFluid::FindRoot_fdf(real_t &root, gsl_function_fdf gsl_func, gsl_root_fdfsolver *s){
    gsl_root_fdfsolver_set (s, &gsl_func, root);
    int status;
    real_t epsrel = 3e-3;
    real_t epsabs = 0;
    len_t max_iter = 30;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        gsl_root_fdfsolver_iterate (s);
        real_t root_prev = root;
        root    = gsl_root_fdfsolver_root (s);
        status = gsl_root_test_delta(root, root_prev, epsabs, epsrel);
        if (status == GSL_SUCCESS)
            break;
    }
}


/**
 * A (crude) method which expands the interval [x_lower, x_upper] 
 * if a root of the provided gsl_function is not in the interval.
 */
void RunawayFluid::FindInterval(real_t *x_lower, real_t *x_upper, gsl_function gsl_func ){
    bool isLoUnderestimate = (gsl_func.function(*x_lower, gsl_func.params ) > 0);
    bool isUpOverestimate = (gsl_func.function(*x_upper, gsl_func.params ) < 0);
    while(!isLoUnderestimate){
        *x_upper = *x_lower;
        *x_lower *= 0.7;
        isLoUnderestimate = (gsl_func.function(*x_lower, gsl_func.params ) > 0);
        isUpOverestimate = true;
    }
    while (!isUpOverestimate){
        *x_lower = *x_upper;
        *x_upper *= 1.4;
        isUpOverestimate = (gsl_func.function(*x_upper, gsl_func.params ) < 0);
    }
}

///////////////////////////////////////////////////////////////
/////////// BEGINNING OF BLOCK WITH METHODS RELATED ///////////
/////////// TO THE CALCULATION OF GROWTH RATES      /////////// 
///////////////////////////////////////////////////////////////

/**
 * Calculates and stores runaway growth rates, such as the avalanche growth. 
 * Uses the matched formula from Hesslow et al. NF 59, 084004 (2019) for
 * the critical runaway momentum, which has been generalized to account for 
 * arbitrary inhomogeneous magnetic fields, see DREAM/doc/notes/theory.
 */
void RunawayFluid::CalculateGrowthRates(){
    real_t *E      = unknowns->GetUnknownData(id_Eterm);
    real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    real_t *n_tot  = unknowns->GetUnknownData(id_ntot); 
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);

    for (len_t ir = 0; ir<this->nr; ir++){
        avalancheGrowthRate[ir] = n_tot[ir] * constPreFactor * criticalREMomentumInvSq[ir];
        real_t pc = criticalREMomentum[ir]; 
        tritiumRate[ir] = evaluateTritiumRate(pc);
        comptonRate[ir] = evaluateComptonRate(pc, compton_photon_flux, gsl_ad_w);
        DComptonRateDpc[ir] = evaluateDComptonRateDpc(pc,compton_photon_flux, gsl_ad_w);

        // Dreicer runaway rate
        bool nnapp = false;
        if (dreicer_nn != nullptr)
            nnapp = dreicer_nn->IsApplicable(T_cold[ir]);  // Is neural network applicable?

        // Neural network
        if (nnapp && (dreicer_mode == OptionConstants::EQTERM_DREICER_MODE_NEURAL_NETWORK
                   || dreicer_mode == OptionConstants::EQTERM_DREICER_MODE_NONE))
            dreicerRunawayRate[ir] = dreicer_nn->RunawayRate(ir, E[ir], n_tot[ir], T_cold[ir]);

        // Connor-Hastie formula
        else {
            real_t Zeff = this->ions->GetZeff(ir);
            dreicerRunawayRate[ir] = dreicer_ConnorHastie->RunawayRate(ir, E[ir], n_cold[ir], Zeff);

            // Emit warning if the Connor-Hastie is the fallback method because
            // we're outside the range of validity of the neural network.
            if (not nnapp && dreicer_mode == OptionConstants::EQTERM_DREICER_MODE_NEURAL_NETWORK)
                DREAM::IO::PrintWarning(
                    DREAM::IO::WARNING_DREICER_NEURAL_NETWORK_INVALID,
                    "Temperature is outside the range of validity for the neural network. "
                    "Falling back to the Connor-Hastie formula instead."
                );
        } 
    }
}

/**
 * Returns the normalized runaway rate due to beta decay of tritium. The net 
 * runaway rate dnRE/dt is obtained after multiplication by n_tritium.
 */
real_t RunawayFluid::evaluateTritiumRate(real_t pc){
    if(isinf(pc))
        return 0;
    real_t gamma_c = sqrt(1+pc*pc);
    real_t gammaMinusOne = pc*pc/(gamma_c+1);
    real_t w = Constants::mc2inEV * gammaMinusOne / tritiumDecayEnergyEV;
    real_t fracAbovePc = 1 - sqrt(w)*w*( 35.0/8.0 - w*21.0/4.0 + w*w*15.0/8.0);
    if(fracAbovePc < 0)
        return 0;

    return log(2) /tritiumHalfLife * fracAbovePc;
}

/**
 * Returns the total cross section for Compton scattering into p>pc due to incident photons 
 * of energy Eg (units of mc and mc2). Eq (29) in Martin-Solis NF 2017.
 */
real_t RunawayFluid::evaluateComptonTotalCrossSectionAtP(real_t Eg, real_t pc){
    real_t gamma_c = sqrt(1+pc*pc);
    real_t x = Eg;
    real_t x2 = x*x;
    real_t x3 = x2*x;
    real_t Wc = pc*pc/(gamma_c+1); // = gamma_c-1
    real_t cc = 1 - 1/Eg * Wc /( Eg - Wc );
    real_t r = 1+x*(1-cc);
    return M_PI * Constants::r0 * Constants::r0 * ( (x2-2*x-2)/x3 * log( (1+2*x)/r ) 
        + 1/(2*x) * ( 1/(r*r) - 1/( (1+2*x)*(1+2*x) ) ) 
        - 1/x3 * ( 1 - x - (1+2*x) / r - x*cc )   );
}

real_t RunawayFluid::evaluateDSigmaComptonDpcAtP(real_t Eg, real_t pc){
    real_t gamma_c = sqrt(1+pc*pc);
    real_t x = Eg;
    real_t x2 = x*x;
    real_t x3 = x2*x;
    real_t Wc = pc*pc/(gamma_c+1); // = gamma_c-1
    real_t cc = 1 - 1/Eg * Wc /( Eg - Wc );
    real_t r = 1+x*(1-cc);
    return M_PI * Constants::r0 * Constants::r0 * ( 
        - (x2-2*x-2)/x3 *  x/(1+2*x)                        // dSigma_compton/d(cosTheta_c)
        +   1/(r*r*r) + 1/x3 * ( (1+2*x)*x / (r*r) + x )  
        ) * (-1/x*(1/(x-Wc)-Wc/x2/((1-Wc/x)*(1-Wc/x))))     // d(cosTheta_c)/dWc       
        * pc/gamma_c;                                       // dWc/dpc                                
}

// Integral of the photon flux spectrum over all Eg (in units of mc2).
const real_t NORMALIZATION_INTEGRATED_COMPTON_SPECTRUM = 5.8844;
/**
 * Returns the photon spectral flux density expected for ITER, Eq (24) in Martin-Solis NF 2017.
 */
real_t RunawayFluid::evaluateComptonPhotonFluxSpectrum(real_t Eg, real_t photonFlux){
    real_t z = (1.2 + log(Eg * Constants::mc2inEV/1e6) ) / 0.8;
    return photonFlux * exp( - exp(-z) - z + 1 ) / NORMALIZATION_INTEGRATED_COMPTON_SPECTRUM;
}

/**
 * Returns the integrand appearing in the evaluation of the total production rate integral (flux density x cross section ) 
 */
struct ComptonParam {real_t pc; real_t photonFlux;};
real_t ComptonIntegrandFunc(real_t Eg, void *par){
    struct ComptonParam *params = (struct ComptonParam *) par;
    
    real_t pc = params->pc;
    real_t photonFlux = params->photonFlux;

    return RunawayFluid::evaluateComptonPhotonFluxSpectrum(Eg, photonFlux) * RunawayFluid::evaluateComptonTotalCrossSectionAtP(Eg,pc);
}

/**
 * Returns the integrand appearing in the evaluation of the derivative w r t pc of
 * the total production rate integral (d/dpc( flux density x cross section )) 
 */
real_t DComptonDpcIntegrandFunc(real_t Eg, void *par){
    struct ComptonParam *params = (struct ComptonParam *) par;
    
    real_t pc = params->pc;
    real_t photonFlux = params->photonFlux;

    return RunawayFluid::evaluateComptonPhotonFluxSpectrum(Eg, photonFlux) * RunawayFluid::evaluateDSigmaComptonDpcAtP(Eg,pc);
}

/**
 * Returns the runaway rate due to Compton scattering on gamma rays. The net runaway rate
 * dnRE/dt is obtained after multiplication by the total electron density n_tot.
 */
real_t RunawayFluid::evaluateComptonRate(real_t pc, real_t photonFlux, gsl_integration_workspace *gsl_ad_w){
    if(isinf(pc))
        return 0;
    real_t gamma_c = sqrt(1+pc*pc);
    real_t gammacMinusOne = pc*pc/(gamma_c+1); // = gamma_c-1
    struct ComptonParam  params= {pc, photonFlux};
    gsl_function ComptonFunc;
    ComptonFunc.function = &(ComptonIntegrandFunc);
    ComptonFunc.params = &params;

    real_t Eg_min = (pc + gammacMinusOne) /2;
    real_t valIntegral;
    // qagiu assumes an infinite upper boundary
    real_t epsrel = 1e-4;
    real_t error;
    gsl_integration_qagiu(&ComptonFunc, Eg_min , 0, epsrel, gsl_ad_w->limit, gsl_ad_w, &valIntegral, &error);
    return valIntegral;
}

/**
 * Returns the derivative of the runaway rate due to Compton scattering on gamma rays w r t pc (factor n_tot NOT included). 
 */
real_t RunawayFluid::evaluateDComptonRateDpc(real_t pc,real_t photonFlux, gsl_integration_workspace *gsl_ad_w){
    if(isinf(pc))
        return 0;
    real_t gamma_c = sqrt(1+pc*pc);
    real_t gammacMinusOne = pc*pc/(gamma_c+1); // = gamma_c-1
    struct ComptonParam  params = {pc, photonFlux};
    gsl_function ComptonFunc;
    ComptonFunc.function = &(DComptonDpcIntegrandFunc);
    ComptonFunc.params = &params;

    real_t Eg_min = (pc + gammacMinusOne) /2;
    real_t valIntegral;
    // qagiu assumes an infinite upper boundary
    real_t epsrel = 1e-4;
    real_t error;
    gsl_integration_qagiu(&ComptonFunc, Eg_min , 0, epsrel, gsl_ad_w->limit, gsl_ad_w, &valIntegral, &error);
    return valIntegral;
}

/**
 * Parameter struct used for the evaluation of pStarFunction.
 */
struct pStarFuncParams {real_t constTerm; len_t ir; RunawayFluid *rf;CollisionQuantity::collqty_settings *collSettingsForPc;};

/**
 * Returns the value of the function whose root (with respect to momentum p) 
 * corresponds to the critical runaway momentum.
 */
real_t RunawayFluid::pStarFunction(real_t p, void *par){
    struct pStarFuncParams *params = (struct pStarFuncParams *) par;
    CollisionQuantity::collqty_settings *collSettingsForPc = params->collSettingsForPc;
    real_t constTerm = params->constTerm;
    real_t ir = params->ir;
    RunawayFluid *rf = params->rf;
    real_t barNuS = rf->evaluateNuSHat(ir,p,collSettingsForPc);
    real_t barNuD = rf->evaluateNuDHat(ir,p,collSettingsForPc);
    return sqrt(sqrt(barNuS*(barNuD+4*barNuS)))/constTerm -  p; 
}

/**
 * Returns the value of the function whose root (with respect to momentum p) 
 * corresponds to the critical runaway momentum.
 */
real_t RunawayFluid::pStarFunctionAlt(real_t p, void *par){
    struct pStarFuncParams *params = (struct pStarFuncParams *) par;
    CollisionQuantity::collqty_settings *collSettingsForPc = params->collSettingsForPc;
    real_t constTerm = params->constTerm;
    real_t ir = params->ir;
    RunawayFluid *rf = params->rf;
    real_t barNuS = rf->evaluateNuSHat(ir,p,collSettingsForPc);
    real_t barNuD = rf->evaluateNuDHat(ir,p,collSettingsForPc);
    return sqrt(sqrt(barNuS*barNuD))/constTerm -  p;
}

/**
 * Calculates pStar with a root finding algorithm for 
 * a given electric field E and radial grid point ir.
 */
real_t RunawayFluid::evaluatePStar(len_t ir, real_t E, gsl_function gsl_func, real_t *nuSHat_COMPSCREEN){
    real_t pStar;
    // Estimate bounds on pStar assuming the limits of complete and no screening. 
    // Note that nuSHat and nuDHat are here independent of p (except via Coulomb logarithm)
    CollisionQuantity::collqty_settings collSetCompScreen = *collSettingsForPc;
    collSetCompScreen.collfreq_type = OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED;
    CollisionQuantity::collqty_settings collSetNoScreen = *collSettingsForPc;
    collSetNoScreen.collfreq_type = OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED;

    *nuSHat_COMPSCREEN = evaluateNuSHat(ir,1,&collSetCompScreen);
    real_t nuDHat_COMPSCREEN = evaluateNuDHat(ir,1,&collSetCompScreen);
    real_t nuSHat_NOSCREEN = evaluateNuSHat(ir,1,&collSetNoScreen);
    real_t nuDHat_NOSCREEN = evaluateNuDHat(ir,1,&collSetNoScreen);

    pc_COMPLETESCREENING[ir] = sqrt(sqrt(*nuSHat_COMPSCREEN*(nuDHat_COMPSCREEN+4**nuSHat_COMPSCREEN))/E);
    pc_NOSCREENING[ir] = sqrt( sqrt(nuSHat_NOSCREEN*(nuDHat_NOSCREEN+4*nuSHat_NOSCREEN)) /E );

    real_t pLo = pc_COMPLETESCREENING[ir];
    real_t pUp = pc_NOSCREENING[ir];
    FindInterval(&pLo,&pUp, gsl_func);
    FindRoot(pLo,pUp, &pStar, gsl_func,fsolve);

    return pStar;
}

/**
 * Calculates and stores the critical runaway momentum. We separately store 1/p^2, since this is the factor
 * entering the avalanche growth rate, and our model will allow it to go negative to capture runaway decay.
 */
void RunawayFluid::CalculateCriticalMomentum(){
    real_t E, constTerm;
    gsl_function gsl_func;
    pStarFuncParams pStar_params;
    real_t pStar;
    real_t nuSHat_COMPSCREEN;
    real_t nuSnuDTerm;
    real_t *E_term = unknowns->GetUnknownData(id_Eterm); 
    for(len_t ir=0; ir<this->nr; ir++){
        /**
         * The normalized electric field E is to be used in the determination of
         * pStar: it is not allowed to be smaller than Eceff in order to behave
         * well in the limit E->0.
         */
        if(fabs(E_term[ir]) > effectiveCriticalField[ir])
            E =  Constants::ec * fabs(E_term[ir]) /(Constants::me * Constants::c);
        else
            E =  Constants::ec * effectiveCriticalField[ir] /(Constants::me * Constants::c);

        real_t EMinusEceff = Constants::ec * (fabs(E_term[ir]) - effectiveCriticalField[ir]) /(Constants::me * Constants::c);

        /**
         * Chooses whether trapping effects are accounted for in growth rates via setting 
         * (could imagine another setting where you go smoothly from one to the other as 
         * t_orbit/t_coll_at_pstar goes from <<1 to >>1)
         */
        real_t effectivePassingFraction = 1;
        if(collSettingsForPc->pstar_mode == OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS)
            effectivePassingFraction = rGrid->GetEffPassFrac(ir);

        constTerm = sqrt(sqrt(E*E * effectivePassingFraction));

        pStar_params = {constTerm,ir,this, collSettingsForPc}; 
        gsl_func.params = &pStar_params;

        if(ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID_HESSLOW){
            gsl_func.function = &(pStarFunctionAlt);
            pStar = evaluatePStar(ir, E, gsl_func, &nuSHat_COMPSCREEN);

            real_t s = pStar*constTerm;
            nuSnuDTerm = s*s*s*s + 4*nuSHat_COMPSCREEN*nuSHat_COMPSCREEN;
        } else {
	        gsl_func.function = &(pStarFunction);
	        pStar = evaluatePStar(ir, E, gsl_func, &nuSHat_COMPSCREEN);

            real_t s = pStar*constTerm;
	        nuSnuDTerm = s*s*s*s;
        }
        
        // Set 1/pc^2 which is to be used in the avalanche growth rate which contains this factor;
        // note that it is allowed to be negative for E<Eceff
        criticalREMomentumInvSq[ir] = EMinusEceff*sqrt(effectivePassingFraction) / sqrt(nuSnuDTerm);

        // also store pc for use in other source functions, but which for E<Eceff is set to inf.
        if (EMinusEceff<=0)
            criticalREMomentum[ir] = std::numeric_limits<real_t>::infinity() ; // should make growth rates zero
        else
            criticalREMomentum[ir] = 1/sqrt(criticalREMomentumInvSq[ir]);
    }
}
    
/**
 *  Returns nuS*p^3/gamma^2, which is constant for ideal plasmas. (only lnL energy dependence)
 */
real_t RunawayFluid::evaluateNuSHat(len_t ir, real_t p, CollisionQuantity::collqty_settings *inSettings){
    return constPreFactor * nuS->evaluateAtP(ir,p,inSettings) / nuS->evaluatePreFactorAtP(p, inSettings->collfreq_mode);
}
/** 
 * Returns nuD*p^3/gamma, which is constant for ideal plasmas. (only lnL energy dependence)
 */
real_t RunawayFluid::evaluateNuDHat(len_t ir, real_t p, CollisionQuantity::collqty_settings *inSettings){
    return constPreFactor * nuD->evaluateAtP(ir,p,inSettings) / nuD->evaluatePreFactorAtP(p, inSettings->collfreq_mode);
}

/**
 * Allocate all stored quantities. 
 */
void RunawayFluid::AllocateQuantities(){
    DeallocateQuantities();

    Ec_free  = new real_t[nr];
    Ec_tot   = new real_t[nr];
    tauEERel = new real_t[nr];
    tauEETh  = new real_t[nr];
    EDreic   = new real_t[nr];

    effectiveCriticalField  = new real_t[nr]; 
    criticalREMomentum      = new real_t[nr];
    criticalREMomentumInvSq = new real_t[nr];
    pc_COMPLETESCREENING    = new real_t[nr];
    pc_NOSCREENING          = new real_t[nr];
    avalancheGrowthRate     = new real_t[nr];
    dreicerRunawayRate      = new real_t[nr];

    tritiumRate = new real_t[nr];
    comptonRate = new real_t[nr];
    DComptonRateDpc = new real_t[nr];

    electricConductivity = new real_t[nr];
}

/**
 * Deallocate all quantities
 */
void RunawayFluid::DeallocateQuantities(){
    if(Ec_free != nullptr){
        delete [] Ec_free;
        delete [] Ec_tot;
        delete [] tauEERel;
        delete [] tauEETh;
        delete [] EDreic;
        delete [] effectiveCriticalField;
        delete [] criticalREMomentum;
        delete [] criticalREMomentumInvSq;
        delete [] pc_COMPLETESCREENING;
        delete [] pc_NOSCREENING;
        delete [] avalancheGrowthRate;
        delete [] dreicerRunawayRate;
        delete [] tritiumRate;
        delete [] comptonRate;
        delete [] DComptonRateDpc;
        delete [] electricConductivity;
    }
}

/**
 * Returns the Sauter electric conductivity of a relativistic plasma.
 *  sigma_Sauter = ( j_||/B ) / (<E*B>/<B^2>)
 * i.e. sigma_Sauter = j_ohm / (E_term / sqrt(<B^2>/Bmin^2) ) 
 */
real_t RunawayFluid::evaluateSauterElectricConductivity(len_t ir, real_t Tcold, real_t Zeff, real_t ncold, bool collisionless){
    return evaluateBraamsElectricConductivity(ir,Tcold,Zeff) * evaluateNeoclassicalConductivityCorrection(ir,Tcold,Zeff,ncold,collisionless);
}

real_t RunawayFluid::evaluateSauterElectricConductivity(len_t ir, bool collisionless){
    return evaluateSauterElectricConductivity(ir, Tcold[ir], ions->GetZeff(ir), ncold[ir], collisionless);
}

/**
 * Returns the Braams-Karney electric conductivity of a relativistic plasma.
 */
real_t RunawayFluid::evaluateBraamsElectricConductivity(len_t ir, real_t Tcold, real_t Zeff){
    if(Zeff<0)
        throw FVM::FVMException("Conductivity: Negative Zeff provided, aborting.");
    const real_t T_SI = Tcold * Constants::ec;

    real_t sigmaBar = gsl_interp2d_eval(gsl_cond, conductivityTmc2, conductivityX, conductivityBraams, 
                T_SI / (Constants::me * Constants::c * Constants::c), 1.0/(1+Zeff), gsl_xacc, gsl_yacc  );
    
    real_t BraamsConductivity = 4*M_PI*Constants::eps0*Constants::eps0 * T_SI*sqrt(T_SI) / 
            (Zeff * sqrt(Constants::me) * Constants::ec * Constants::ec * lnLambdaEE->GetLnLambdaT(ir) ) * sigmaBar;
    return BraamsConductivity;
}
real_t RunawayFluid::evaluateBraamsElectricConductivity(len_t ir){
    return evaluateBraamsElectricConductivity(ir, Tcold[ir], ions->GetZeff(ir));
}
/**
 * Returns the correction to the Spitzer conductivity, valid in all collisionality regimes,
 * taken from O Sauter, C Angioni and Y R Lin-Liu, Phys Plasmas 6, 2834 (1999).
 */
real_t RunawayFluid::evaluateNeoclassicalConductivityCorrection(len_t ir, real_t Tcold, real_t Zeff, real_t ncold, bool collisionLess){
    real_t ft = 1 - rGrid->GetEffPassFrac(ir);
    
    real_t X = ft;
    const real_t R0 = rGrid->GetR0();
    if(isinf(R0))
        X = 0;
    else if(!collisionLess){
        // qR0 is the safety factor multiplied by R0
        const real_t *jtot = unknowns->GetUnknownData(id_jtot);
        real_t mu0Ip = Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir,rGrid,jtot);
        const real_t qR0 = rGrid->SafetyFactorNormalized(ir,mu0Ip);
        real_t TkeV = Tcold/1000;
        real_t eps = rGrid->GetR(ir)/R0;
        real_t nuEStar = 0.012*(ncold/1e20)*Zeff * qR0/(eps*sqrt(eps) * TkeV*TkeV);

        X /= 1 + (0.55-0.1*ft)*sqrt(nuEStar) + 0.45*(1-ft)*nuEStar/(Zeff*sqrt(Zeff)) ;
    }
    return 1 - (1+0.36/Zeff)*X + X*X/Zeff * (0.59-0.23*X);
}

real_t RunawayFluid::evaluateNeoclassicalConductivityCorrection(len_t ir, bool collisionLess){
    return evaluateNeoclassicalConductivityCorrection(ir, Tcold[ir], ions->GetZeff(ir), ncold[ir], collisionLess);
}

/**
 * Returns the partial derivative of the conductivity with respect to unknown derivId,
 * choosing conductivity formula based on cond_mode setting
 */
real_t RunawayFluid::evaluatePartialContributionConductivity(len_t ir, len_t derivId, len_t n){
    switch(cond_mode) {
        case OptionConstants::CONDUCTIVITY_MODE_BRAAMS: 
            return evaluatePartialContributionBraamsConductivity(ir, derivId, n);
            break;
        case OptionConstants::CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS: 
            return evaluatePartialContributionSauterConductivity(ir, derivId, n, true);
            break;
        case OptionConstants::CONDUCTIVITY_MODE_SAUTER_COLLISIONAL: 
            return evaluatePartialContributionSauterConductivity(ir, derivId, n, false);
            break;
        default:
            return std::numeric_limits<real_t>::infinity();
            break; // throw exception?            
    }
}

/**
 * Calculation of the partial derivative of Sauter conductivity
 */  
real_t RunawayFluid::evaluatePartialContributionSauterConductivity(len_t ir, len_t derivId, len_t n, bool collisionless) {
    real_t eps = std::numeric_limits<real_t>::epsilon();
    real_t Zeff = ions->GetZeff(ir);
    if(derivId==id_Tcold){
        real_t h = Tcold[ir]*sqrt(eps);
        return ( evaluateSauterElectricConductivity(ir,Tcold[ir]+h,Zeff,ncold[ir],collisionless)
               - evaluateSauterElectricConductivity(ir,Tcold[ir]-h,Zeff,ncold[ir],collisionless) ) / (2*h);
    } else if (derivId==id_ncold) {
        real_t h = ncold[ir]*sqrt(eps);
        return ( evaluateSauterElectricConductivity(ir,Tcold[ir],Zeff,ncold[ir]+h,collisionless)
               - evaluateSauterElectricConductivity(ir,Tcold[ir],Zeff,ncold[ir]-h,collisionless) ) / (2*h);
    } else if (derivId==id_ni){
        // using dZeff/dni = Z0^2/nfree - Z0*<Z0^2>/nfree^2
        len_t iz,Z0;
        ions->GetIonIndices(n,iz,Z0);
        real_t nfree = ions->GetFreeElectronDensityFromQuasiNeutrality(ir);
        if(nfree==0)
            return 0;
        real_t nZ0Z0 = ions->GetNZ0Z0(ir);
        real_t h = 1e-6*Zeff;
        return Z0/nfree * (Z0 - nZ0Z0/nfree) * 
            ( evaluateSauterElectricConductivity(ir,Tcold[ir],Zeff+h,ncold[ir],collisionless)
            - evaluateSauterElectricConductivity(ir,Tcold[ir],Zeff-h,ncold[ir],collisionless) ) / (2*h);
    } else 
        return 0;
}

/**
 * Calculation of the partial derivative of Sauter conductivity
 */  
real_t RunawayFluid::evaluatePartialContributionBraamsConductivity(len_t ir, len_t derivId, len_t n) {
    real_t eps = std::numeric_limits<real_t>::epsilon();
    real_t Zeff = ions->GetZeff(ir);
    if(derivId==id_Tcold){
        real_t h = Tcold[ir]*sqrt(eps);
        return ( evaluateBraamsElectricConductivity(ir,Tcold[ir]+h,Zeff)
               - evaluateBraamsElectricConductivity(ir,Tcold[ir]-h,Zeff) ) / (2*h);
    } else if (derivId==id_ni){
        // using dZeff/dni = Z0^2/nfree - Z0*<Z0^2>/nfree^2
        len_t iz,Z0;
        ions->GetIonIndices(n,iz,Z0);
        real_t nfree = ions->GetFreeElectronDensityFromQuasiNeutrality(ir);
        if(nfree==0)
            return 0;
        real_t nZ0Z0 = ions->GetNZ0Z0(ir);
        real_t h = 1e-6*Zeff;
        return Z0/nfree * (Z0 - nZ0Z0/nfree) * 
            ( evaluateBraamsElectricConductivity(ir,Tcold[ir],Zeff+h)
            - evaluateBraamsElectricConductivity(ir,Tcold[ir],Zeff-h) ) / (2*h);
    } else 
        return 0;
}

/**
 * Calculation of the partial derivative of the avalanche growth rate 
 * with respect to unknown quantities. So far approximate expression
 * assuming the E-field dependence is captured via the (E-Eceff) coefficient
 * and density via Eceff ~ n_tot.
 */
void RunawayFluid::evaluatePartialContributionAvalancheGrowthRate(real_t *dGamma, len_t derivId) {
    if( !( (derivId==id_Eterm) || (derivId==id_ntot) ) ){
        for(len_t ir = 0; ir<nr; ir++)
            dGamma[ir] = 0;
    }else{
        // set dGamma to d(Gamma)/d(E_term)
        for(len_t ir=0; ir<nr; ir++)
            dGamma[ir] = avalancheGrowthRate[ir] / ( fabs(Eterm[ir]) - effectiveCriticalField[ir] );

        // if derivative w.r.t. n_tot, multiply by d(E-Eceff)/dntot = -dEceff/dntot ~ -Eceff/ntot
        if(derivId==id_ntot)
            for(len_t ir=0; ir<nr; ir++)
                dGamma[ir] *= - effectiveCriticalField[ir] / ntot[ir];
        // else multiply by sign of E
        else if(derivId==id_Eterm)
            for(len_t ir=0;ir<nr;ir++){
                real_t sgnE = (Eterm[ir]>0) - (Eterm[ir]<0);
                dGamma[ir] *= sgnE;
            }
    }
}

/**
 * Calculation of the partial derivative of the compton scattering growth rate 
 * with respect to unknown quantities, assuming pc~sqrt(ntot/(E-Eceff)). Note 
 * also that although Eg_min depends on pc, the cross section is zero at Eg_min.
 */
void RunawayFluid::evaluatePartialContributionComptonGrowthRate(real_t *dGamma, len_t derivId) {

    if(derivId==id_Eterm){
        // set dGamma to d(Gamma)/d(E_term)
        for(len_t ir=0; ir<nr; ir++){
            if(isinf(criticalREMomentum[ir]))
                dGamma[ir]=0;
            else {
                real_t sgnE = (Eterm[ir]>0) - (Eterm[ir]<0);
                dGamma[ir] = -1/2* DComptonRateDpc[ir] * criticalREMomentum[ir] * sgnE/( fabs(Eterm[ir]) - effectiveCriticalField[ir] ) ;
            }
        }
    } else if (derivId==id_ntot){
        // set dGamma to d(Gamma)/d(ntot)-gamma_compton/ntot
        // NOTE! This only includes the effect of ntot on pc, and not the explicit derivative!
        // The explicit derivative will be added automatically since the compton source is implemented
        // as a DiagonalComplexTerm
        for(len_t ir=0; ir<nr; ir++){
            if(isinf(criticalREMomentum[ir]))
                dGamma[ir]=0;
            else
                dGamma[ir] = 1/2* DComptonRateDpc[ir] * criticalREMomentum[ir]/ntot[ir] ;
        }
    } else
        for(len_t ir = 0; ir<nr; ir++)
            dGamma[ir] = 0;
}

/**
 * Printing timing information for this object.
 */
void RunawayFluid::PrintTimings() {
    this->timeKeeper->PrintTimings(true, timerTot);
}

/**
 * Save timing information to the given SFile object.
 *
 * sf:   SFile object to save timing info to.
 * path: Path in SFile object to save timing info to.
 */
void RunawayFluid::SaveTimings(SFile *sf, const std::string& path) {
    this->timeKeeper->SaveTimings(sf, path);
}
