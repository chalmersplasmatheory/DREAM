/**
 * Implementation of collision-rate calculator that calculates
 * various collision, ionisation, recombination, growth etc rates and quantities.  
 * It takes an UnknownQuantityHandler, extracts needed parameters, loads atomic  
 * physics data and calculates a bunch of collison-related quantities.
 * Also allows manual specification of plasma parameters (or even collision frequencies etc).
*/


/** 
 * EXAMPLE: DREAM simulation workflow
 * 
 * // initialize
 * CollisionQuantityHandler *cqh = new CollisionQuantityHandler(struct collqtyhand_settings*);
 * CollQty->SetGrid(grid);
 * CollQty->SetUnknowns(unknowns);
 * 
 * // each time plasma parameters have changed, update collision rates etc:
 * CollQty->Rebuild(); // calculates collision frequencies, ionisation rates and derived quantities (growth rates etc)
 */



#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>




using namespace DREAM;

const len_t CollisionQuantityHandler::ionSizeAj_len = 55; 
const real_t CollisionQuantityHandler::ionSizeAj_data[ionSizeAj_len] = { 0.631757734322417, 0.449864664424796, 0.580073385681175, 0.417413282378673, 0.244965367639212, 0.213757911761448, 0.523908484242040, 0.432318176055981, 0.347483799585738, 0.256926098516580, 0.153148466772533, 0.140508604177553, 0.492749302776189, 0.419791849305259, 0.353418389488286, 0.288707775999513, 0.215438905215275, 0.129010899184783, 0.119987816515379, 0.403855887938967, 0.366602498048607, 0.329462647492495, 0.293062618368335, 0.259424839110224, 0.226161504309134, 0.190841656429844, 0.144834685411878, 0.087561370494245, 0.083302176729104, 0.351554934261205, 0.328774241757188, 0.305994557639981, 0.283122417984972, 0.260975850956140, 0.238925715853581, 0.216494264086975, 0.194295316086760, 0.171699132959493, 0.161221485564969, 0.150642403738712, 0.139526182041846, 0.128059339783537, 0.115255069413773, 0.099875435538094, 0.077085983503479, 0.047108093547224, 0.045962185039177, 0.235824746357894, 0.230045911002090, 0.224217341261303, 0.215062179624586, 0.118920957451653, 0.091511805821898, 0.067255603181663, 0.045824624741631 };
const real_t CollisionQuantityHandler::ionSizeAj_Zs[ionSizeAj_len] = { 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 54, 54, 54, 74, 74, 74, 74, 74 };
const real_t CollisionQuantityHandler::ionSizeAj_Z0s[ionSizeAj_len] = { 0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 1, 2, 3, 0, 30, 40, 50, 60 };

const len_t CollisionQuantityHandler::meanExcI_len = 39;
const real_t CollisionQuantityHandler::meanExcI_data[meanExcI_len] = { 8.3523e-05, 1.1718e-04, 6.4775e-05, 2.1155e-04, 2.6243e-04, 1.2896e-04, 1.8121e-04, 
        2.6380e-04, 4.1918e-04, 9.5147e-04, 0.0011, 2.6849e-04, 3.2329e-04, 3.8532e-04, 4.6027e-04, 5.5342e-04, 
        6.9002e-04, 9.2955e-04, 0.0014, 0.0028, 0.0029, 3.6888e-04, 4.2935e-04, 4.9667e-04, 5.7417e-04, 6.6360e-04, 
        7.7202e-04, 9.0685e-04, 0.0011, 0.0014, 0.0016, 0.0017, 0.0019, 0.0022, 0.0027, 0.0035, 0.0049, 0.0092, 0.0095};
const real_t CollisionQuantityHandler::meanExcI_Zs[meanExcI_len] = { 2, 2, 3, 3, 3, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 10, 
        10, 10, 10, 10, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18};
const real_t CollisionQuantityHandler::meanExcI_Z0s[meanExcI_len] = { 0, 1, 0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 
        4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};

/** 
 * Constructor
 */ 
CollisionQuantityHandler::CollisionQuantityHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  enum OptionConstants::momentumgrid_type mgtype,  struct collqtyhand_settings *cqset){
    ionHandler = ih;
    grid       = g;
    unknowns   = u;
    settings   = cqset;
    gridtype   = mgtype;
}

/**
 * Destructor.
 */
CollisionQuantityHandler::~CollisionQuantityHandler(){
    DeallocateCollisionFrequencies();
    DeallocateIonisationRates();
    DeallocateLnLambdas();
    //DeallocateIonSpecies();
    DeallocateDerivedQuantities();
    DeallocateHiGi();
    DeallocateGSL();

}

/**
 * Calculates and stores all collision quantities from an UnknownQuantityHandler. 
 */
void CollisionQuantityHandler::Rebuild() {

    len_t id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->n_cold   = unknowns->GetUnknownData(id_ncold);

    len_t id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->T_cold   = unknowns->GetUnknownData(id_Tcold);
   
    this->nZ = ionHandler->GetNZ();
    this->nzs = ionHandler->GetNzs();
//    this->ionDensity     = ionHandler->GetDensityMat();
    this->ZAtomicNumber  = ionHandler->GetZs();
//    this->Z0ChargeNumber = ionHandler->GetZ0List();
    
    /*
    n_tot = new real_t[n];
    for (len_t ir=0; ir < n; ir++){
        for (len_t iz = 0; iz < nZ; iz++){
            n_tot[ir] += this->ionDensity[ir][iz] * ZAtomicNumber[iz];
        }
    }
    */

    /*
    bool gridHasChanged;
    bool temperatureHasChanged = unknowns->HasChanged(id_Tcold);
    bool densitiesHaveChanged  = unknowns->HasChanged(GetUnknownID(OptionConstants::UQTY_ION_SPECIES));

    if (gridHasChanged)
        CalculateGridQuantities();


    if (temperatureHasChanged){
        if (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
            InitializeGSLWorkspace();
    }
    if (temperatureHasChanged || densitiesHaveChanged || gridHasChanged){
        CalculateCoulombLogarithms(); 
//        CalculateHcoldGcold();
        CalculateHiGiFuncs();                    // hi, gi, hcold, gcold
        CalculateCollisionFrequenciesFromHiGi(); // nu_s, nu_D and nu_||
    }
    */

    

    /**
     * The following three methods calculate and store all related functions, 
     * which may reappear in later calculations (for example in building the 
     * Newton-method Jacobian matrix). By instead running 
     * CalculateCollisionFrequencies(); we would only store nu_s, nu_D and nu_||. 
     */ 
    CalculateCoulombLogarithms();            // all lnLs. Rename to CalculateThermalQuantities and include hcold, gcold? Then hi, gi, the heavy parts, only need to be rebuilt if grid changes 
    CalculateHiGiFuncs();                    // hi, gi, hcold, gcold
    CalculateCollisionFrequenciesFromHiGi(); // nu_s, nu_D and nu_||

    CalculateIonisationRates();

    CalculateDerivedQuantities();
}


/**
 * Initializes a GSL workspace for each radius (used for relativistic test particle operator evaluation),
 * using a T_cold-dependent fixed quadrature. 
 */
void CollisionQuantityHandler::InitializeGSLWorkspace(){
 /** 
  * (consider using a single regular dynamic quadrature instead as the integral is somewhat tricky, 
  * since in the limit p/mc -> 0 the integral is sharply peaked at p_min -- goes as int 1/sqrt(x) dx,0,inf --
  * and may be challenging to resolve using a fixed point quadrature)
  */
    DeallocateGSL();
    gsl_w = new gsl_integration_fixed_workspace*[n];
    const real_t lowerLim = 0; // integrate from 0 to inf
    const gsl_integration_fixed_type *T = gsl_integration_fixed_laguerre;
    const len_t Npoints = 20; // play around with this number -- may require larger, or even sufficient with lower
    const real_t alpha = 0.0;
    real_t b;
    real_t Theta;
    for (len_t ir = 0; ir<n; ir++){
        Theta = T_cold[ir]/Constants::mc2inEV;
        b = 1/Theta;
        gsl_w[ir] = gsl_integration_fixed_alloc(T, Npoints, lowerLim, b, alpha, 0.0);
    }
}

/**
 * Calculates n_cold contribution to nu_s
 */
real_t CollisionQuantityHandler::evaluateHColdAtP(len_t i, real_t p) {    
    // Depending on setting, set nu_s to superthermal or full formula (with maxwellian)
    if (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL)
        return evaluateLnLambdaEEAtP(i,p) * constPreFactor * (1+p*p)/(p*p*p);
    else if (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        // nu_s = lnLee * constPreFactor * M / p^3;
        // M = (gamma^2 * Psi1 - Theta*Psi0 + (Theta*gamma-1)*p*exp(- (gamma-1)/Theta) )/[ exp(1/Theta)K_2(1/Theta) ];
        // Psi0 = int_0^p exp( -(sqrt(1+s^2)-1)/Theta) / sqrt(1+s^2) ds;
        // Psi1 = int_0^p exp( -(sqrt(1+s^2)-1)/Theta) ds;
        real_t gamma = sqrt(1+p*p);
        real_t Theta = T_cold[i] / Constants::mc2inEV;
        real_t M = 0;
        M += gamma*gamma* evaluatePsi1(Theta,p) - Theta * evaluatePsi0(Theta,p);
        M +=  (Theta*gamma - 1) * p * exp( -(gamma-1)/Theta );
        M /= evaluateExp1OverThetaK(Theta,2.0);
        return evaluateLnLambdaEEAtP(i,p) * constPreFactor * M / (p*p*p);
    } else
        throw NotImplementedException("Chosen collfreq_mode setting not yet supported.");
}




/**
 *  Calculates ion contribution to nu_s
 */
real_t CollisionQuantityHandler::evaluateHiAtP(len_t i, real_t p, len_t Z, len_t Z0) {    
    
    if (settings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED)
        return evaluateBetheHiAtP(p,Z,Z0);
    else if (settings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED)
        return (Z-Z0) * evaluateHColdAtP(i, p);
    else 
        return 0;

}

/**
 * Calculates n_cold contribution to nu_D
 */
real_t CollisionQuantityHandler::evaluateGColdAtP(len_t i, real_t p) {
    // Depending on setting, set nu_D to superthermal or full formula (with maxwellian)
    real_t p2 = p*p;
    real_t gamma = sqrt(1+p2);
    real_t Theta;
    real_t M = 0;
    if (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL)
        return evaluateLnLambdaEEAtP(i,p) * constPreFactor * gamma/(p2*p);
    else if (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        Theta = T_cold[i] / Constants::mc2inEV;
        M += (p2*gamma*gamma + Theta*Theta)*evaluatePsi0(i,p);
        M += Theta*(2*p2*p2 - 1)*evaluatePsi1(i,p);
        M += gamma*Theta * ( 1 + Theta*(2*p2-1)*p*exp( -(gamma-1)/Theta ) );
        M /= evaluateExp1OverThetaK(Theta,2.0);
        return evaluateLnLambdaEEAtP(i,p) * constPreFactor * M  / (gamma * p2*p2*p);
    } else
        throw NotImplementedException("Chosen collfreq_mode setting not yet supported.");


}

/**
 * Calculates ion contribution to nu_D
 */
real_t CollisionQuantityHandler::evaluateGiAtP(len_t i, real_t p, len_t Z, len_t Z0) {
    real_t g_i;
    real_t lnL = evaluateLnLambdaEIAtP(i,p);
    real_t constBit =  constPreFactor * sqrt(1+p*p)/(p*p*p);
    // the completely screened contribution
    g_i = Z0*Z0 * lnL * constBit;
    
    if (settings->collfreq_type == OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED)
        g_i += (Z*Z-Z0*Z0) * lnL * constBit;
    else if (settings->collfreq_type == OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED){
        g_i += evaluateKirillovGiAtP(p, Z, Z0);
    }
    
    return g_i;
}



/** 
 * Evaluates partially screened ion contribution to nu_D using Kirillov's model
 */
real_t CollisionQuantityHandler::evaluateKirillovGiAtP(real_t p, len_t Z, len_t Z0){
    // Using DFT-obtained (where available) or analytical a_j.
    real_t aj = GetIonEffectiveSizeAj(Z,Z0);
    real_t x = pow( p*aj, 3/2 );
    return constPreFactor * sqrt(1+p*p)/(p*p*p) * (2/3) * (
            (Z*Z - Z0*Z0)*log( 1+x ) -  (Z-Z0)*(Z-Z0) * x/(1+x) );
}


/**
 * Evaluates partially screened ion contribution to nu_s using Bethe's formula
 */
real_t CollisionQuantityHandler::evaluateBetheHiAtP(real_t p, len_t Z, len_t Z0){
    real_t gamma = sqrt(p*p+1);
    real_t h = p*sqrt(gamma-1) / GetMeanExcitationEnergy(Z,Z0);
    real_t k = 5;
    // if (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL)
    //    return constPreFactor * (Z-Z0) * gamma*gamma/(p*p*p) * ( log( h ) - p*p/(gamma*gamma) );
    //else if (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) 
        // This one should probably always be used? 
        return  constPreFactor * (Z-Z0) * gamma*gamma/(p*p*p) * ( log(1+ pow(h,k)) / k  - p*p/(gamma*gamma) ) ;
    
    //else 
    //    return -1; // no such setting implemented

}

real_t CollisionQuantityHandler::GetIonEffectiveSizeAj(len_t Z, len_t Z0){
    // Fetch DFT-calculated value from table if it exists:
    for (len_t n=0; n<ionSizeAj_len; n++)
        if( Z==ionSizeAj_Zs[n] && (Z0==ionSizeAj_Z0s[n]) )
            return ionSizeAj_data[n];

    // If DFT-data is missing, use Kirillov's model:
    return 2/Constants::alpha * pow(9*M_PI,1/3) / 4 * pow(Z-Z0,2/3) / Z;
}

real_t CollisionQuantityHandler::GetMeanExcitationEnergy(len_t Z, len_t Z0){
    // Fetch value from table if it exists:
    for (len_t n=0; n<meanExcI_len; n++)
        if( Z==meanExcI_Zs[n] && (Z0==meanExcI_Z0s[n]) )
            return meanExcI_data[n];

    // if can't find in the table, return something large so that the contribution 
    // to nu_s becomes zero. Or send an error? 
    return __DBL_MAX__; 

}


real_t CollisionQuantityHandler::psi0Integrand(real_t x, void *params){
    real_t gamma = *(real_t *) params;
    return 1/sqrt( (x+gamma)*(x+gamma)-1 );
} 
real_t CollisionQuantityHandler::psi1Integrand(real_t x, void *params){
    real_t gamma = *(real_t *) params;
    return (x+gamma)/sqrt((x+gamma)*(x+gamma)-1); // integrated with weight w(x) = exp(-(x-gamma)/Theta) 
} 
/** 
 * Evaluates integral appearing in relativistic test-particle operator
 * Psi0 = int_0^p exp( -(sqrt(1+s^2)-1)/Theta) / sqrt(1+s^2) ds;
 */
real_t CollisionQuantityHandler::evaluatePsi0(len_t ir, real_t p) {
    real_t gamma = sqrt(1+p*p);
    
    
    gsl_function F;
    F.function = &(CollisionQuantityHandler::psi0Integrand); 
    F.params = &gamma;
    real_t psi0int; 
    gsl_integration_fixed(&F, &psi0int, gsl_w[ir]);

    
    real_t Theta = T_cold[ir] / Constants::mc2inEV;
    return evaluateExp1OverThetaK(Theta,0) - exp( -(gamma-1)/Theta ) * psi0int;

}
real_t CollisionQuantityHandler::evaluatePsi1(len_t ir, real_t p) {
    
    real_t gamma = sqrt(1+p*p);
    gsl_function F;
    F.function = &(CollisionQuantityHandler::psi1Integrand); 
    F.params = &gamma;
    real_t psi1int; 
    gsl_integration_fixed(&F, &psi1int, gsl_w[ir]);

    real_t Theta = T_cold[ir] / Constants::mc2inEV;
    return evaluateExp1OverThetaK(Theta,1) - exp( -(gamma-1)/Theta ) * psi1int;
    

}


real_t CollisionQuantityHandler::evaluateExp1OverThetaK(real_t Theta, real_t n) {
    real_t ThetaThreshold = 0.002;
    /**
     * Since cyl_bessel_k ~ exp(-1/Theta), for small Theta you get precision issues.
     * Instead using asymptotic expansion for bessel_k for small Theta.
     */
    if (Theta > ThetaThreshold)
        return exp(1/Theta)*std::cyl_bessel_k(n,1/Theta);
    else {
//        return sqrt(M_PI*Theta/2)*(1 + 15*Theta/8 + 105*Theta*Theta/128 - 945*Theta*Theta*Theta/3072);
        real_t n2 = n*n;
        return sqrt(M_PI*Theta/2)*(1 + (4*n2-1)/8 * Theta + (4*n2-1)*(4*n2-9)*Theta*Theta/128 + (4*n2-1)*(4*n2-9)*(4*n2-25)*Theta*Theta*Theta/3072);
    }
}


/**
 * Calculates and stores nu_s and nu_D
 */
void CollisionQuantityHandler::CalculateCollisionFrequenciesFromHiGi(){
    DeallocateCollisionFrequencies();
    real_t 
        **nu_s    = nullptr, //new real_t*[n], 
        **nu_D    = nullptr, //new real_t*[n], 
        **nu_par  = nullptr, //new real_t*[n],
        **nu_D2   = new real_t*[n],
        **nu_s1   = new real_t*[n], 
        **nu_s2   = new real_t*[n], 
        **nu_D1   = new real_t*[n],
        **nu_par1 = new real_t*[n],
        **nu_par2 = new real_t*[n];
    real_t /*p,*/ p_f1, p_f2;
    bool collfreqmodeFull = (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);

    len_t ind;
    for (len_t ir = 0; ir < n; ir++) {
        
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        /* For now safe not to calculate or store anything on distribution grid
        nu_s[ir]   = new real_t[np1*np2];
        nu_D[ir]   = new real_t[np1*np2];
        nu_par[ir] = new real_t[np1*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                p = mg->GetP(i,j);
                ind = j*np1+i;

                nu_s[ir][ind] = n_cold[ir]*HCold[ir][ind];
                nu_D[ir][ind] = n_cold[ir]*GCold[ir][ind];
                for (len_t iz = 0; iz<nZ[ir]; iz++) {
                    nu_s[ir][ind] += ionDensity[ir][iz] * HiFunc[ir][ind][iz];
                    nu_D[ir][ind] += ionDensity[ir][iz] * GiFunc[ir][ind][iz];
                }
                if (collfreqmodeFull)
                    nu_par[ir][ind] = nu_s[ir][ind] * (T_cold[ir]/Constants::mc2inEV)*sqrt(1+p*p);
                else 
                    // not sure if there is a middle ground where superthermal nu_par could be interesting.. 
                    // I don't think that it would work (inf particle flux at p=0 for maxwellian)
                    nu_par[ir][ind] = 0;
            }
        }
        */

        // Terms on flux grid 1
        nu_s1[ir]   = new real_t[(np1+1)*np2];
        nu_D1[ir]   = new real_t[(np1+1)*np2];
        nu_par1[ir] = new real_t[(np1+1)*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                p_f1 = mg->GetP_f1(i,j);
                ind = j*(np1+1)+i;

                nu_s1[ir][ind] = n_cold[ir]*HCold_f1[ir][ind];
                nu_D1[ir][ind] = n_cold[ir]*GCold_f1[ir][ind];
                for (len_t iz = 0; iz<nZ; iz++) {
                    for (len_t Z0 = 0; Z0<ZAtomicNumber[iz]+1; Z0++){
                        nu_s1[ir][ind] += ionHandler->GetIonDensity(ir,iz,Z0) * HiFunc_f1[ir][ind][iz];
                        nu_D1[ir][ind] += ionHandler->GetIonDensity(ir,iz,Z0) * GiFunc_f1[ir][ind][iz];
                    }
                }
                if (collfreqmodeFull)
                    nu_par1[ir][ind] = nu_s1[ir][ind] * (T_cold[ir]/Constants::mc2inEV)*sqrt(1+p_f1*p_f1);
                else 
                    nu_par1[ir][ind] = 0;
            }
        }

        // Terms on flux grid 2
        nu_s2[ir]   = new real_t[np1*(np2+1)];
        nu_D2[ir]   = new real_t[np1*(np2+1)];
        nu_par2[ir] = new real_t[np1*(np2+1)];
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                p_f2 = mg->GetP_f2(i,j);
                ind = j*np1+i;

                nu_s2[ir][ind] = n_cold[ir]*HCold_f2[ir][ind];
                nu_D2[ir][ind] = n_cold[ir]*GCold_f2[ir][ind];
                for (len_t iz = 0; iz<nZ; iz++) {
                    for (len_t Z0 = 0; Z0<ZAtomicNumber[iz]+1; Z0++){
                        nu_s2[ir][ind] += ionHandler->GetIonDensity(ir,iz,Z0) * HiFunc_f2[ir][ind][iz];
                        nu_D2[ir][ind] += ionHandler->GetIonDensity(ir,iz,Z0) * GiFunc_f2[ir][ind][iz];
                    }
                }
                if (collfreqmodeFull)
                    nu_par2[ir][ind] = nu_s2[ir][ind] * (T_cold[ir]/Constants::mc2inEV)*sqrt(1+p_f2*p_f2);
                else 
                    nu_par2[ir][ind] = 0;
            }
        }
    }
    this->collisionFrequencyNuS      = nu_s;
    this->collisionFrequencyNuS_f1   = nu_s1;
    this->collisionFrequencyNuS_f2   = nu_s2;
    this->collisionFrequencyNuD      = nu_D;
    this->collisionFrequencyNuD_f1   = nu_D1;
    this->collisionFrequencyNuD_f2   = nu_D2;
    this->collisionFrequencyNuPar    = nu_par;
    this->collisionFrequencyNuPar_f1 = nu_par1;
    this->collisionFrequencyNuPar_f2 = nu_par2;
    
}


/**
 * Calculates and stores h_i and g_i -- the partial contributions to nu_s and nu_D from ions i --  on n x nz x (np1 x np2).
 */
void CollisionQuantityHandler::CalculateHiGiFuncs(){
    DeallocateHiGi();
    real_t 
        ***hi       = nullptr, //new real_t**[n],
         **gCold    = nullptr, //new real_t*[n],
        ***gi       = nullptr, //new real_t**[n],
         **hCold    = nullptr, //new real_t*[n],
        ***hi_f1    = new real_t**[n],
        ***hi_f2    = new real_t**[n],
        ***gi_f1    = new real_t**[n],
        ***gi_f2    = new real_t**[n],
         **hCold_f1 = new real_t*[n],
         **hCold_f2 = new real_t*[n],
         **gCold_f1 = new real_t*[n],
         **gCold_f2 = new real_t*[n];

    for (len_t ir=0; ir<n; ir++){
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        len_t np1 = mg->GetNp1();
        len_t np2 = mg->GetNp2();
        real_t /*p,*/ p_f1, p_f2;
        
        /* For now safe not to calculate or store anything on distribution grid
        hi[ir]    = new real_t*[np1*np2];
        gi[ir]    = new real_t*[np1*np2];
        hCold[ir] = new real_t[np1*np2];
        gCold[ir] = new real_t[np1*np2];

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                p = mg->GetP(i,j);    
                
                hCold[ir][j*np1+i] = evaluateHColdAtP(ir,p);
                gCold[ir][j*np1+i] = evaluateGColdAtP(ir,p);
                
                hi[ir][j*np1+i] = new real_t[nZ[ir]];
                gi[ir][j*np1+i] = new real_t[nZ[ir]];
                for (len_t iz=0; iz<nZ[ir]; iz++){
                    hi[ir][j*np1+i][iz] = evaluateHiAtP(ir,p,ZAtomicNumber[ir][iz],Z0ChargeNumber[ir][iz]);
                    gi[ir][j*np1+i][iz] = evaluateGiAtP(ir,p,ZAtomicNumber[ir][iz],Z0ChargeNumber[ir][iz]);
                }
            }
        }
        */

        hi_f1[ir]    = new real_t*[(np1+1)*np2];
        gi_f1[ir]    = new real_t*[(np1+1)*np2];
        hCold_f1[ir] = new real_t[(np1+1)*np2];
        gCold_f1[ir] = new real_t[(np1+1)*np2];
        
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                p_f1 = mg->GetP_f1(i,j);
                
                hCold_f1[ir][j*(np1+1)+i] = evaluateHColdAtP(ir,p_f1);
                gCold_f1[ir][j*(np1+1)+i] = evaluateGColdAtP(ir,p_f1);
                
                hi_f1[ir][j*(np1+1)+i] = new real_t[nzs];
                gi_f1[ir][j*(np1+1)+i] = new real_t[nzs];
                for (len_t iz=0; iz<nZ; iz++){
                    for (len_t Z0 = 0; Z0<ZAtomicNumber[iz]+1; Z0++){
                        hi_f1[ir][j*(np1+1)+i][ionHandler->GetIndex(iz,Z0)] = evaluateHiAtP(ir,p_f1,ZAtomicNumber[iz],Z0);
                        gi_f1[ir][j*(np1+1)+i][ionHandler->GetIndex(iz,Z0)] = evaluateGiAtP(ir,p_f1,ZAtomicNumber[iz],Z0);    
                    }
                }
            }
        }
        
        gi_f2[ir]    = new real_t*[np1*(np2+1)];
        hi_f2[ir]    = new real_t*[np1*(np2+1)];
        gCold_f2[ir] = new real_t[np1*(np2+1)];
        hCold_f2[ir] = new real_t[np1*(np2+1)];
        
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                p_f2 = mg->GetP_f2(i,j);

                hCold_f2[ir][j*np1+i] = evaluateHColdAtP(ir,p_f2);
                gCold_f2[ir][j*np1+i] = evaluateGColdAtP(ir,p_f2);
                
                hi_f2[ir][j*np1+i] = new real_t[nzs];
                gi_f2[ir][j*np1+i] = new real_t[nzs];
                for (len_t iz=0; iz<nZ; iz++){
                    for (len_t Z0 = 0; Z0<ZAtomicNumber[iz]+1; Z0++){
                        hi_f2[ir][j*np1+i][ionHandler->GetIndex(iz,Z0)] = evaluateHiAtP(ir,p_f2,ZAtomicNumber[iz],Z0);
                        gi_f2[ir][j*np1+i][ionHandler->GetIndex(iz,Z0)] = evaluateGiAtP(ir,p_f2,ZAtomicNumber[iz],Z0);
                    }
                }
            }
        }
    }
    this->HiFunc    = hi;
    this->HiFunc_f1 = hi_f1;
    this->HiFunc_f2 = hi_f2;
    this->GiFunc    = gi;
    this->GiFunc_f1 = gi_f1;
    this->GiFunc_f2 = gi_f2;
    this->HCold     = hCold;
    this->HCold_f1  = hCold_f1;
    this->HCold_f2  = hCold_f2;
    this->GCold     = gCold;
    this->GCold_f1  = gCold_f1;
    this->GCold_f2  = gCold_f2;
}



// Bonus function not used by the DREAM simulation workflow
void CollisionQuantityHandler::CalculateCoulombLogarithms(){ 
    DeallocateLnLambdas();
    //if (ionDensity==nullptr){
        // error?
    //}
    real_t /* p, gamma,*/ p_f1, p_f2,  gamma_f1, gamma_f2;
    
    real_t 
        **lnLee  = nullptr, //new real_t*[n], 
        **lnLei  = nullptr, //new real_t*[n], 
        **lnLee1 = new real_t*[n], 
        **lnLee2 = new real_t*[n], 
        **lnLei1 = new real_t*[n], 
        **lnLei2 = new real_t*[n],
         *lnLc   = new real_t[n], 
         *lnLTe  = new real_t[n];


    for (len_t ir = 0; ir < n; ir++) {
        
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        lnLc[ir]  = 14.6 + 0.5*log( T_cold[ir]/(n_cold[ir]/1e20) );
        lnLTe[ir] = 14.9 + 0.5*log( (T_cold[ir]/1e3)*(T_cold[ir]/1e3)/(n_cold[ir]/1e20) );

        /* For now safe not to calculate or store anything on distribution grid
        lnLee[ir] = new real_t[np1*np2];
        lnLei[ir] = new real_t[np1*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                p     = mg->GetP(i,j);
                gamma = sqrt(1+p*p);
                lnLee[ir][j*np1+i] = lnLc[ir] + log( sqrt(gamma-1) );
                lnLei[ir][j*np1+i] = lnLc[ir] + log( sqrt(2)*p );
            }
        }
        */

        lnLee1[ir] = new real_t[(np1+1)*np2];
        lnLei1[ir] = new real_t[(np1+1)*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                p_f1 = mg->GetP_f1(i,j);
                gamma_f1 = sqrt(1+p_f1*p_f1);
                lnLee1[ir][j*(np1+1)+i] = lnLc[ir] + log( sqrt(gamma_f1-1) );
                lnLei1[ir][j*(np1+1)+i] = lnLc[ir] + log( sqrt(2)*p_f1 );
            }
        }

        lnLee2[ir] = new real_t[np1*(np2+1)];
        lnLei2[ir] = new real_t[np1*(np2+1)];
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                p_f2 = mg->GetP_f2(i,j);
                gamma_f2 = sqrt(1+p_f2*p_f2);
                lnLee2[ir][j*np1+i] = lnLc[ir] + log( sqrt(gamma_f2-1) );
                lnLei2[ir][j*np1+i] = lnLc[ir] + log( sqrt(2)*p_f2 );
            }
        }
    }

    this->lnLambda_c     = lnLc;
    this->lnLambda_Te    = lnLTe;
    this->lnLambda_ee    = lnLee;
    this->lnLambda_ee_f1 = lnLee1;
    this->lnLambda_ee_f2 = lnLee2;
    this->lnLambda_ei    = lnLei;
    this->lnLambda_ei_f1 = lnLei1;
    this->lnLambda_ei_f2 = lnLei2;
    
}



// Loads ADAS coefficients and uses ion species and atomic parmeters data to calculate
//   various ionisation and recombination rates
void CollisionQuantityHandler::CalculateIonisationRates(){
    DeallocateIonisationRates();

  
    // SetIonisationRates(Icold, Ikin, IRE,
    //                    RR, CEZP, CEHP){
    
}




// Uses collision frequencies and ion species to calculate
// critical fields and avalanche growth rates
void CollisionQuantityHandler::CalculateDerivedQuantities(){
    Ec_free = new real_t[n];
    Ec_tot  = new real_t[n];

    for (len_t ir=0; ir<n; ir++){
        Ec_free[ir] = lnLambda_c[ir] * n_cold[ir] * constPreFactor * Constants::me * Constants::c / Constants::ec;
        Ec_tot[ir]  = lnLambda_c[ir] * n_tot[ir]  * constPreFactor * Constants::me * Constants::c / Constants::ec;
    }
    CalculateEffectiveCriticalField();
    CalculatePStar();

    CalculateGrowthRates();

    DeallocateDerivedQuantities();
    // SetDerivedQuantities(Ec, Ectot, ED, 
    //                        Gamma_avalanche, Eceff);
    
}







void CollisionQuantityHandler::DeallocateCollisionFrequencies(){
    if (this->collisionFrequencyNuD_f2 == nullptr)
        return;


    for (len_t i = 0; i < n; i++) {
        delete [] this->collisionFrequencyNuS[i];
        delete [] this->collisionFrequencyNuS_f1[i];
        delete [] this->collisionFrequencyNuS_f2[i];
        delete [] this->collisionFrequencyNuD[i];
        delete [] this->collisionFrequencyNuD_f1[i];
        delete [] this->collisionFrequencyNuD_f2[i];
        delete [] this->collisionFrequencyNuPar[i];
        delete [] this->collisionFrequencyNuPar_f1[i];
        delete [] this->collisionFrequencyNuPar_f2[i];
        
    }

    delete [] this->collisionFrequencyNuS;
    delete [] this->collisionFrequencyNuS_f1;
    delete [] this->collisionFrequencyNuS_f2;
    delete [] this->collisionFrequencyNuD;
    delete [] this->collisionFrequencyNuD_f1;
    delete [] this->collisionFrequencyNuD_f2;
    delete [] this->collisionFrequencyNuPar;
    delete [] this->collisionFrequencyNuPar_f1;
    delete [] this->collisionFrequencyNuPar_f2;
    
}


void CollisionQuantityHandler::DeallocateHiGi(){
    if (this->GiFunc_f2 == nullptr)
        return;

    for (len_t ir = 0; ir < n; ir++) {
        delete [] this->HCold[ir];
        delete [] this->GCold[ir];
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        len_t np1 = mg->GetNp1();
        len_t np2 = mg->GetNp2();
    
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                delete [] this->HiFunc[ir][j*np1+i];
                delete [] this->GiFunc[ir][j*np1+i];
            }
        }
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                delete [] this->HiFunc_f1[ir][j*(np1+1)+i];
                delete [] this->GiFunc_f1[ir][j*(np1+1)+i];
            }
        }
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                delete [] this->HiFunc_f2[ir][j*np1+i];
                delete [] this->GiFunc_f2[ir][j*np1+i];
            }
        }
        delete [] this->HiFunc[ir];
        delete [] this->HiFunc_f1[ir];
        delete [] this->HiFunc_f2[ir];
        delete [] this->GiFunc[ir];
        delete [] this->GiFunc_f1[ir];
        delete [] this->GiFunc_f2[ir];

    }

    delete [] this->HiFunc;
    delete [] this->HiFunc_f1;
    delete [] this->HiFunc_f2;
    delete [] this->GiFunc;
    delete [] this->GiFunc_f1;
    delete [] this->GiFunc_f2;
    delete [] this->HCold;
    delete [] this->GCold;
}






void CollisionQuantityHandler::DeallocateLnLambdas(){
    if (this->lnLambda_c == nullptr)
        return;

    delete [] this->lnLambda_c;
    delete [] this->lnLambda_Te;
    

    for (len_t i = 0; i < n; i++) {
        delete [] this->lnLambda_ee[i];
        delete [] this->lnLambda_ee_f1[i];
        delete [] this->lnLambda_ee_f2[i];
        delete [] this->lnLambda_ei[i];
        delete [] this->lnLambda_ei_f1[i];
        delete [] this->lnLambda_ei_f2[i];
    }

    delete [] this->lnLambda_ee;
    delete [] this->lnLambda_ee_f1;
    delete [] this->lnLambda_ee_f2;
    delete [] this->lnLambda_ei;
    delete [] this->lnLambda_ei_f1;
    delete [] this->lnLambda_ei_f2;
}




void CollisionQuantityHandler::DeallocateGSL(){
    if (this->gsl_w == nullptr)
        return;

    for (len_t ir=0; ir<n; ir++)
        gsl_integration_fixed_free(gsl_w[ir]);
}


void CollisionQuantityHandler::DeallocateIonisationRates(){
    if (this->ionisationRateCold == nullptr)
        return;

    for (len_t ir = 0; ir<n; ir++){
        delete [] this->ionisationRateCold[ir];
        delete [] this->ionisationRateHot[ir];
        delete [] this->ionisationRateREFluid[ir];
        delete [] this->recombinationRateRadiative[ir];
        delete [] this->chargeExchangeZP[ir];
    }
    delete [] this->ionisationRateCold;
    delete [] this->ionisationRateHot;
    delete [] this->ionisationRateREFluid;
    delete [] this->recombinationRateRadiative;
    delete [] this->chargeExchangeZP;
    delete [] this->chargeExchangeHP;
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
real_t UFrictionTermIntegrand(real_t xi0, void *par){
    struct UFuncParams *params = (struct UFuncParams *) par;
    len_t ir = params->ir;
    real_t A = params->A;
    FVM::RadialGrid *rGrid = params->rGrid;
    std::function<real_t(real_t,real_t,real_t)> FrictionTermFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return xi0/sqrt(1-BOverBmin*(1-xi0*xi0))*BOverBmin;};
    return rGrid->CalculateFluxSurfaceAverage(ir,false, FrictionTermFunc)*exp(-A*(1-xi0));
}

real_t USynchrotronTermIntegrand(real_t xi0, void *par){
    struct UFuncParams *params = (struct UFuncParams *) par;
    len_t ir = params->ir;
    real_t A = params->A;
    FVM::RadialGrid *rGrid = params->rGrid;
    std::function<real_t(real_t,real_t,real_t)> SynchrotronTermFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return (1-xi0*xi0)*sqrt(1-BOverBmin*(1-xi0*xi0))*BOverBmin*BOverBmin*BOverBmin;};
    return rGrid->CalculateFluxSurfaceAverage(ir,false, SynchrotronTermFunc)*exp(-A*(1-xi0));
}


// Evaluates the effective momentum flow U accounting for electric field, collisional friction and radiation reaction 
real_t CollisionQuantityHandler::evaluateUAtP(len_t ir,real_t p, real_t Eterm,gsl_integration_workspace *gsl_ad_w){
    FVM::RadialGrid *rGrid =  grid->GetRadialGrid();
    const real_t Bmin = rGrid->GetBmin(ir);
    const real_t Bmax = rGrid->GetBmax(ir);
    const real_t B2avg = rGrid->GetFSA_B2(ir);
    
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avg)/Bmin; 
    real_t xiT = sqrt(1-Bmin/Bmax);
    real_t A = 2*E/(p*evaluateNuDAtP(ir,p));
    
    real_t Econtrib = E/(A*A) *( A-1 * exp(-A*(1-xiT))*(A*xiT -1) );

    real_t FrictionTerm = p*evaluateNuSAtP(ir,p);
    /* Uncomment when we can support bremsstrahlung losses
    if(!(settings->bremsstrahlung_mode==OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT))
        FrictionTerm += evaluateBremsStoppingForceAtP(ir,p) / (Constants::me * Constants::c);
    */
    UFuncParams FuncParams = {ir, A, rGrid};
    gsl_function UIntegrandFunc;

    UIntegrandFunc.function = &(UFrictionTermIntegrand);
    UIntegrandFunc.params = &FuncParams;
    real_t frictionIntegral;
    gsl_integration_qags(&UIntegrandFunc, xiT,1.0,0,1e-4,1000,gsl_ad_w, &frictionIntegral, nullptr);
    real_t FrictionContrib = -FrictionTerm * frictionIntegral;

    real_t SynchrotronTerm = Constants::ec * Constants::ec * Constants::ec * Constants::ec * Bmin * Bmin
                            / ( 6 * M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me
                                * Constants::c * Constants::c * Constants::c * sqrt(1+p*p));
    UIntegrandFunc.function = &(USynchrotronTermIntegrand);
    real_t synchrotronIntegral;
    gsl_integration_qags(&UIntegrandFunc, xiT,1.0,0,1e-4,1000,gsl_ad_w, &synchrotronIntegral, nullptr);
    real_t SynchrotronContrib = -SynchrotronTerm * synchrotronIntegral;
    
   return Econtrib + FrictionContrib + SynchrotronContrib;
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
    gsl_min_fminimizer *s;

    // Requiring that the solution lies between p=1 and p=500... anything else would be very unusual
    real_t p_ex_guess = 10.0;
    real_t p_ex_lo = 1.0, p_ex_up = 500.0;
    gsl_function F;

    CollisionQuantityHandler::UExtremumParams params = {ir,Eterm,gsl_ad_w,collQtyHand};
    F.function = &(UExtremumFunc);
    F.params = &params;

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc(T);
    gsl_min_fminimizer_set(s, &F, p_ex_guess, p_ex_lo, p_ex_up);


    int status;
    real_t rel_error = 1e-3;
    len_t max_iter = 10;
    for (len_t iteration = 0; iteration < max_iter; iteration++ ){
        status  = gsl_min_fminimizer_iterate(s);
        *p_ex   = gsl_min_fminimizer_x_minimum(s);
        p_ex_lo = gsl_min_fminimizer_x_lower(s);
        p_ex_up = gsl_min_fminimizer_x_upper(s);
        status  = gsl_root_test_interval(p_ex_lo, p_ex_up, 0, rel_error);

        if (status == GSL_SUCCESS){
            gsl_min_fminimizer_free(s);
            break;
        }
    }

    // Return function value -U(p_ex)
    return gsl_min_fminimizer_f_minimum(s);
}

// For GSL function: returns minimum of -U at given Eterm -- This could be compactified by writing FindUExtremumAtE on this form from the beginning
real_t ECritFunc(real_t E, void *par){
    struct CollisionQuantityHandler::UExtremumParams *params = (struct CollisionQuantityHandler::UExtremumParams *) par;
    len_t ir = params->ir;
    gsl_integration_workspace *gsl_w = params->gsl_w;
    CollisionQuantityHandler *collQtyHand = params->collQtyHand;

    return FindUExtremumAtE(ir, E, nullptr, gsl_w,collQtyHand);    
}


/** 
 * The function evaluates the effective critical field, defined (in doc/notes/theory.pdf) as the minimum electric field
 * (or rather, Eterm, effectively a loop voltage) above which U(p)=0 has real roots in p, where U is the pitch
 * averaged momentum flux in the limit (E-Eceff)<<E. It can be a bit hard to penetrate due to all GSL function thrown around,
 * but essentially it performs a root finding algorithm to solve the problem U_max(Eceff) = 0. U_max is in turn obtained using
 * a minimization algorithm in -U(p; E), where finally U(p) is obtained using gsl integration of various flux surface averages.
 * It performs up to 100 (max_iter*max_iter) function evaluations of nu_s and nu_D for each radial index ir, and in each such evaluation also evaluates 
 * multiple flux surface averages (inside a gsl adaptive quadrature, so tens?).
 */
void CollisionQuantityHandler::CalculateEffectiveCriticalField(){
    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);
    
    //real_t Eceff_guess;

    real_t ELo, EUp;
    UExtremumParams params;
    gsl_function UExtremumFunc;
    for (len_t ir=0; ir<n; ir++){
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
    criticalREMomentum = new real_t[n];

    len_t id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    real_t *E_term = unknowns->GetUnknownData(id_Eterm);

    real_t E, constTerm;
    gsl_function gsl_func;
    pStarFuncParams pStar_params;
    real_t pLo, pUp, pStar;
    for(len_t ir=0; ir<n; ir++){
        if(E_term[ir] > effectiveCriticalField[ir])
            E =  Constants::ec * E_term[ir] /(Constants::me * Constants::c);
        else
            E =  Constants::ec * effectiveCriticalField[ir] /(Constants::me * Constants::c);

        constTerm = sqrt(sqrt(E*E * grid->GetRadialGrid()->GetEffPassFrac(ir)));

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
        criticalREMomentum[ir] =  sqrt(sqrt( (evaluateBarNuSNuDAtP(ir,pStar) + 4*evaluateNuSAtP(ir,pStar)*pStar*pStar*pStar/(1+pStar*pStar))  
                                                / (E*E * grid->GetRadialGrid()->GetEffPassFrac(ir)) ));
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

// Takes a p interval [x_lower, x_upper] and iterates at most max_iter=10 times (or to a relative error of rel_error=0.001)
// to find an estimate for p_Star 
void CollisionQuantityHandler::FindPStarRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func){
    const gsl_root_fsolver_type *GSL_rootsolver_type = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (GSL_rootsolver_type);
    gsl_root_fsolver_set (s, &gsl_func, x_lower, x_upper); 

    int status;
    real_t rel_error = 1e-3;
    len_t max_iter = 10;
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
 * Calculates the runaway rate due to beta decay of tritium. We need to implement a setting for tritiumFraction,
 * since the current ion structure cannot distinguish isotopes -- and it is probably only in the case of tritium
 * that this would be interesting. We could probably have a special input parameter to be tritium fraction (i.e.)
 * fraction of hydrogenic density that is tritium. Elsewhere in the code we can assume the T to behave like D and  
 * that the fraction is constant in radius (? for simplicity at least, but not necessarily).
 */
real_t CollisionQuantityHandler::evaluateTritiumRate(len_t ir){
    /*real_t tritiumFraction=0; 
    real_t nH = 0;
    for(len_t iz=0; iz<nZ; iz++){
        if( (ZAtomicNumber[iz]==1) )
            nH += ionDensity[ir][iz];
    }
    real_t n_tritium = tritiumFraction * nH;
    */
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
    struct ComptonParam  params= {criticalREMomentum[ir], this};
    gsl_function ComptonFunc;
    ComptonFunc.function = &(ComptonIntegrandFunc);
    ComptonFunc.params = &params;

    real_t gamma_c = sqrt(1+criticalREMomentum[ir]*criticalREMomentum[ir]);
    real_t Eg_min = (criticalREMomentum[ir] + gamma_c - 1) /2;
    real_t valIntegral;
    // qagiu assumes an infinite upper boundary
    real_t epsrel = 1e-4;
    gsl_integration_qagiu(&ComptonFunc, Eg_min , 0, epsrel, 1000, gsl_ad_w, &valIntegral, nullptr);
    return n_tot[ir]*valIntegral;
}

// Evaluates and stores the avalanche, tritium and compton growth rates
void CollisionQuantityHandler::CalculateGrowthRates(){
    len_t id_nRE = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
    real_t *nRE = unknowns->GetUnknownData(id_nRE);



    gsl_integration_workspace *gsl_ad_w = gsl_integration_workspace_alloc(1000);


    real_t gamma_crit;
    avalancheRate = new real_t[n];
    tritiumRate   = new real_t[n];
    comptonRate   = new real_t[n];
    for (len_t ir = 0; ir<n; ir++){
        // we still haven't implemented the relativistic corrections in criticalREmomentum, 
        // but let's keep it like this for now in case we do in the future.
        gamma_crit = sqrt( 1 + criticalREMomentum[ir]*criticalREMomentum[ir] );
        avalancheRate[ir] = 0.5 * nRE[ir] * constPreFactor / (gamma_crit-1) ;
        tritiumRate[ir] = evaluateTritiumRate(ir);
        comptonRate[ir] = evaluateComptonRate(ir,gsl_ad_w);
    }
}



//real_t CollisionQuantityHandler::ionSizeAj_data[ionSizeAj_Ntab] = { 8.3523e-05, 1.1718e-04, 6.4775e-05, 2.1155e-04, 2.6243e-04, 1.2896e-04, 1.8121e-04, 2.6380e-04, 4.1918e-04, 9.5147e-04, 0.0011, 2.6849e-04, 3.2329e-04, 3.8532e-04, 4.6027e-04, 5.5342e-04, 6.9002e-04, 9.2955e-04, 0.0014, 0.0028, 0.0029, 3.6888e-04, 4.2935e-04, 4.9667e-04, 5.7417e-04, 6.6360e-04, 7.7202e-04, 9.0685e-04, 0.0011, 0.0014, 0.0016, 0.0017, 0.0019, 0.0022, 0.0027, 0.0035, 0.0049, 0.0092, 0.0095};


/****************************************************
 * Methods not currently used in the DREAM workflow *
 ****************************************************/

/**
 * Calculates and stores nu_s and nu_D
 */
void CollisionQuantityHandler::CalculateCollisionFrequencies(){
    DeallocateCollisionFrequencies();
    real_t 
        **nu_s    = new real_t*[n], 
        **nu_D    = new real_t*[n], 
        **nu_D2   = new real_t*[n],
        **nu_s1   = new real_t*[n], 
        **nu_s2   = new real_t*[n], 
        **nu_D1   = new real_t*[n],
        **nu_par  = new real_t*[n],
        **nu_par1 = new real_t*[n],
        **nu_par2 = new real_t*[n];
    real_t p, p_f1, p_f2;
    bool collfreqmodeFull = (settings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);


    for (len_t ir = 0; ir < n; ir++) {
        
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        // Terms on distribution grid
        nu_s[ir]   = new real_t[np1*np2];
        nu_D[ir]   = new real_t[np1*np2];
        nu_par[ir] = new real_t[np1*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                p = mg->GetP(i,j);
                nu_s[ir][j*np1+i] = evaluateNuSAtP(ir,p);
                nu_D[ir][j*np1+i] = evaluateNuDAtP(ir,p);
                if (collfreqmodeFull)
                    nu_par[ir][j*np1+i] = nu_s[ir][j*np1+i] * (T_cold[ir]/Constants::mc2inEV)*sqrt(1+p*p);
                else // not sure if there is a middle ground where superthermal nu_par could be interesting.. I don't think that it would work (inf particle flux at p=0 for maxwellian)
                    nu_par[ir][j*np1+i] = 0;
            }
        }

        // Terms on flux grid 1
        nu_s1[ir]   = new real_t[(np1+1)*np2];
        nu_D1[ir]   = new real_t[(np1+1)*np2];
        nu_par1[ir] = new real_t[(np1+1)*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                p_f1 = mg->GetP_f1(i,j); 
                nu_s1[ir][j*(np1+1)+i]  = evaluateNuSAtP(ir,p_f1);
                nu_D1[ir][j*(np1+1)+i]  = evaluateNuDAtP(ir,p_f1);
                if (collfreqmodeFull)
                    nu_par1[ir][j*(np1+1)+i] = nu_s[ir][j*(np1+1)+i] * (T_cold[ir]/Constants::mc2inEV)*sqrt(1+p_f1*p_f1);
                else // not sure if there is a middle ground where superthermal nu_par could be interesting.. I don't think that it would work (inf particle flux at p=0 for maxwellian)
                    nu_par1[ir][j*(np1+1)+i] = 0;
            }
        }

        // Terms on flux grid 2
        nu_s2[ir]   = new real_t[np1*(np2+1)];
        nu_D2[ir]   = new real_t[np1*(np2+1)];
        nu_par2[ir] = new real_t[np1*(np2+1)];
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                p_f2 = mg->GetP_f2(i,j);                  
                nu_s2[ir][j*np1+i]  = evaluateNuSAtP(ir,p_f2);
                nu_D2[ir][j*np1+i]  = evaluateNuDAtP(ir,p_f2);
                if (collfreqmodeFull)
                    nu_par2[ir][j*np1+i] = nu_s[ir][j*np1+i] *  (T_cold[ir]/Constants::mc2inEV)*sqrt(1+p_f2*p_f2);
                else // not sure if there is a middle ground where superthermal nu_par could be interesting.. I don't think that it would work (inf particle flux at p=0 for maxwellian)
                    nu_par2[ir][j*np1+i] = 0;
            }
        }
    }
    this->collisionFrequencyNuS      = nu_s;
    this->collisionFrequencyNuS_f1   = nu_s1;
    this->collisionFrequencyNuS_f2   = nu_s2;
    this->collisionFrequencyNuD      = nu_D;
    this->collisionFrequencyNuD_f1   = nu_D1;
    this->collisionFrequencyNuD_f2   = nu_D2;
    this->collisionFrequencyNuPar    = nu_par;
    this->collisionFrequencyNuPar_f1 = nu_par1;
    this->collisionFrequencyNuPar_f2 = nu_par2;
    
}

/*
void CollisionQuantityHandler::DeallocateIonSpecies(){
    if (n_cold == nullptr)
        return;

    delete [] n_cold; 
}

void CollisionQuantityHandler::SetIonSpecies(real_t **dens, len_t *Z, len_t *Z0, real_t *T){
    DeallocateIonSpecies();
    this->ionDensity     = dens;
    this->ZAtomicNumber  = Z;
    //this->Z0ChargeNumber = Z0;
    this->T_cold = T;

    real_t *n_free = new real_t[n];
    for (len_t i   = 0; i<n; i++){
        for (len_t iz = 0; iz<nZ; iz++){
            n_free[i]  += Z0[iz]*dens[i][iz];
        }
    }
    this->n_cold = n_free;
}    
*/






/**
 * Calculates nu_s
 */
real_t CollisionQuantityHandler::evaluateNuSAtP(len_t i, real_t p){
    real_t ns = n_cold[i]*evaluateHColdAtP(i, p);
    
    // sum the contributions from all ion species
    for (len_t iZ = 0; iZ<nZ; iZ++){
        for (len_t Z0=0; Z0<ZAtomicNumber[iZ]+1; Z0++){
            ns += ionHandler->GetIonDensity(i,iZ,Z0)
                * evaluateHiAtP(i, p, ZAtomicNumber[iZ],Z0);
        }
    }

    return ns;
                
}

/**
 * Calculates nu_D
 */
real_t CollisionQuantityHandler::evaluateNuDAtP(len_t i, real_t p){
    real_t nD = n_cold[i]*evaluateGColdAtP(i, p);

    // sum the contributions from all ion species
    for (len_t iZ = 0; iZ<nZ; iZ++){
        for (len_t Z0=0; Z0<ZAtomicNumber[iZ]+1; Z0++){
            nD += ionHandler->GetIonDensity(i,iZ,Z0)
                * evaluateGiAtP(i, p, ZAtomicNumber[iZ],Z0);
        }
    }
                
    return nD;
}


/**
 * Evaluates relativistic constant lnLambda
 */
real_t CollisionQuantityHandler::evaluateLnLambdaC(len_t i) {
    return 14.6 + 0.5*log( T_cold[i]/(n_cold[i]/1e20) );
}

/**
 * Evaluates energy dependent e-e lnLambda
 */
real_t CollisionQuantityHandler::evaluateLnLambdaEEAtP(len_t i, real_t p) {
    real_t gamma = sqrt(p*p+1);

    if (settings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT)
        return evaluateLnLambdaC(i);
    else if (settings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT)
        return evaluateLnLambdaC(i) + log( sqrt(gamma-1) );
    else 
        throw NotImplementedException("Chosen lnL_type setting not yet supported.");
}

/**
 * Evaluates energy dependent e-i lnLambda
 */
real_t CollisionQuantityHandler::evaluateLnLambdaEIAtP(len_t i, real_t p) {
    if (settings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_CONSTANT)
        return evaluateLnLambdaC(i);
    else if (settings->lnL_type==OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT)
        return evaluateLnLambdaC(i) + log( sqrt(2)*p );
    else 
        throw NotImplementedException("Chosen lnL_type setting not yet supported.");
}

