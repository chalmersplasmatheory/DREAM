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

#include <cmath>





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
CollisionQuantityHandler::CollisionQuantityHandler(struct collqtyhand_settings *cq){
    if (cq == nullptr)
        this->settings = new struct collqtyhand_settings;
    else
        this->settings = cq;

    
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

    len_t id_ions  = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    // should use the unknownquantityhandler to do the following 8 lines 
    real_t *ionDensEqSys = unknowns->GetUnknownData(id_ions);
    real_t **ionDensReshaped  = new real_t*[n];

    for (len_t ir=0; ir < n; ir++){
        ionDensReshaped[ir] = new real_t[nZ[ir]];
        for (len_t iz = 0; iz < nZ[ir]; iz++){
            ionDensReshaped[ir][iz] = ionDensEqSys[iz*n+ir]; // the densities are probably to be enumerated this way?
        }
    }
    this->ionDensity = ionDensReshaped;
    //" this->ZAtomicNumber  = unknowns->GetUnknownData(id_ions)->ZAtomicNumber ";
    //" this->Z0ChargeNumber = unknowns->GetUnknownData(id_ions)->Z0ChargeNumber ";
    
    InitializeGSLWorkspace();

    /**
     * The following three methods calculate and store all related functions, 
     * which may reappear in later calculations (for example in building the 
     * Newton-method Jacobian matrix). By instead running 
     * CalculateCollisionFrequencies(); we would only store nu_s, nu_D and nu_||. 
     */ 
    CalculateCoulombLogarithms();            // all lnLs 
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
        return -1;
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
        return -1; // error: no such setting supported

}

/**
 * Calculates ion contribution to nu_D
 */
real_t CollisionQuantityHandler::evaluateGiAtP(len_t i, real_t p, len_t Z, len_t Z0) {
    real_t g_i;

    // the completely screened contribution
    g_i = Z0*Z0 * evaluateLnLambdaEIAtP(i,p) * constPreFactor * sqrt(1+p*p)/(p*p*p);
    
    if (settings->collfreq_type == OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED)
        g_i += (Z*Z-Z0*Z0) * evaluateLnLambdaEIAtP(i,p) * constPreFactor * sqrt(1+p*p)/(p*p*p);
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

    // todo: find/write a quadrature function quad
    // return quad([](real_t s){return exp( -(sqrt(1+s*s)-1)/Theta); }, 0, p, Npoints );

    
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
                for (len_t iz = 0; iz<nZ[ir]; iz++) {
                    nu_s1[ir][ind] += ionDensity[ir][iz] * HiFunc_f1[ir][ind][iz];
                    nu_D1[ir][ind] += ionDensity[ir][iz] * GiFunc_f1[ir][ind][iz];
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
                for (len_t iz = 0; iz<nZ[ir]; iz++) {
                    nu_s2[ir][ind] += ionDensity[ir][iz] * HiFunc_f2[ir][ind][iz];
                    nu_D2[ir][ind] += ionDensity[ir][iz] * GiFunc_f2[ir][ind][iz];
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
                
                hi_f1[ir][j*(np1+1)+i] = new real_t[nZ[ir]];
                gi_f1[ir][j*(np1+1)+i] = new real_t[nZ[ir]];
                for (len_t iz=0; iz<nZ[ir]; iz++){
                    hi_f1[ir][j*(np1+1)+i][iz] = evaluateHiAtP(ir,p_f1,ZAtomicNumber[ir][iz],Z0ChargeNumber[ir][iz]);
                    gi_f1[ir][j*(np1+1)+i][iz] = evaluateGiAtP(ir,p_f1,ZAtomicNumber[ir][iz],Z0ChargeNumber[ir][iz]);    
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
                
                hi_f2[ir][j*np1+i] = new real_t[nZ[ir]];
                gi_f2[ir][j*np1+i] = new real_t[nZ[ir]];
                for (len_t iz=0; iz<nZ[ir]; iz++){    
                    hi_f2[ir][iz][j*np1+i] = evaluateHiAtP(ir,p_f2,ZAtomicNumber[ir][iz],Z0ChargeNumber[ir][iz]);
                    gi_f2[ir][iz][j*np1+i] = evaluateGiAtP(ir,p_f2,ZAtomicNumber[ir][iz],Z0ChargeNumber[ir][iz]);
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
    if (ionDensity==nullptr){
        // error?
    }
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












//real_t CollisionQuantityHandler::ionSizeAj_data[ionSizeAj_Ntab] = { 8.3523e-05, 1.1718e-04, 6.4775e-05, 2.1155e-04, 2.6243e-04, 1.2896e-04, 1.8121e-04, 2.6380e-04, 4.1918e-04, 9.5147e-04, 0.0011, 2.6849e-04, 3.2329e-04, 3.8532e-04, 4.6027e-04, 5.5342e-04, 6.9002e-04, 9.2955e-04, 0.0014, 0.0028, 0.0029, 3.6888e-04, 4.2935e-04, 4.9667e-04, 5.7417e-04, 6.6360e-04, 7.7202e-04, 9.0685e-04, 0.0011, 0.0014, 0.0016, 0.0017, 0.0019, 0.0022, 0.0027, 0.0035, 0.0049, 0.0092, 0.0095};


/** 
 * Methods not currently used in the DREAM workflow
 */

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


void CollisionQuantityHandler::DeallocateIonSpecies(){
    if (n_cold == nullptr)
        return;

    delete [] n_cold; 
}

void CollisionQuantityHandler::SetIonSpecies(real_t **dens, len_t **Z, len_t **Z0, real_t *T){
    DeallocateIonSpecies();
    this->ionDensity     = dens;
    this->ZAtomicNumber  = Z;
    this->Z0ChargeNumber = Z0;
    this->T_cold = T;

    real_t *n_free = new real_t[n];
    for (len_t i   = 0; i<n; i++){
        for (len_t iz = 0; iz<nZ[i]; iz++){
            n_free[i]  += Z0[i][iz]*dens[i][iz];
        }
    }
    this->n_cold = n_free;
}    







/**
 * Calculates nu_s
 */
real_t CollisionQuantityHandler::evaluateNuSAtP(len_t i, real_t p){
    real_t ns = n_cold[i]*evaluateHColdAtP(i, p);
    
    // sum the contributions from all ion species
    for (len_t iZ = 0; iZ<nZ[i]; iZ++)
        ns += ionDensity[i][iZ]
            * evaluateHiAtP(i, p, ZAtomicNumber[i][iZ],Z0ChargeNumber[i][iZ]);

    return ns;
                
}

/**
 * Calculates nu_D
 */
real_t CollisionQuantityHandler::evaluateNuDAtP(len_t i, real_t p){
    real_t nD = n_cold[i]*evaluateGColdAtP(i, p);

    // sum the contributions from all ion species
    for (len_t iZ = 0; iZ<nZ[i]; iZ++)
                    nD += ionDensity[i][iZ]
                        * evaluateGiAtP(i, p, ZAtomicNumber[i][iZ],Z0ChargeNumber[i][iZ]);
                
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
        return -1; //no such setting implemented     
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
        return -1; //no such setting implemented 
}

