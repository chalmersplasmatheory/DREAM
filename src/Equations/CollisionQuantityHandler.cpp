/**
 * Implementation of collision-rate calculator that calculates
 * various collision, ionisation, recombination, growth etc rates and quantities.  
 * It takes an EquationSystem, extracts needed parameters, loads atomic  
 * physics data and calculates a bunch of collison-related quantities.
 * Also allows manual specification of plasma parameters (or even collision frequencies etc).
*/


/** 
 * EXAMPLE: DREAM simulation workflow
 * 
 * // initialize
 * CollisionQuantityHandler *CollQty(collqty_settings);
 * CollQty->SetGrid(grid);
 * CollQty->SetEqSys(EqSys);
 * 
 * // each time plasma parameters have changed, update collision rates etc:
 * CollQty->RebuildFromEqSys(); // calculates collision frequencies, ionisation rates and derived quantities (growth rates etc)
 */



#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/EquationSystem.hpp"

#include <cmath>





using namespace DREAM;

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
    DeallocateIonSpecies();
    DeallocateDerivedQuantities();
}


void CollisionQuantityHandler::RebuildFromEqSys() {

    len_t id_ncold = eqSys->GetUnknownID(SimulationGenerator::UQTY_N_COLD);
    this->n_cold = eqSys->GetUnknownData(id_ncold);

    len_t id_Tcold = eqSys->GetUnknownID(SimulationGenerator::UQTY_T_COLD);
    this->T_cold = eqSys->GetUnknownData(id_Tcold);

    len_t id_ions = eqSys->GetUnknownID(SimulationGenerator::UQTY_ION_SPECIES);
    //this->ionDensity     = eqSys->GetUnknownData(id_ions)->ionDensity;
    //this->ZAtomicNumber  = eqSys->GetUnknownData(id_ions)->ZAtomicNumber;
    //this->Z0ChargeNumber = eqSys->GetUnknownData(id_ions)->Z0ChargeNumber;
    
    /**
     * if (grid or Zs or Z0s changed)
     *     LoadAtomicData();
     */

    CalculateIonisationRates();
    CalculateCollisionFrequencies();
    CalculateDerivedQuantities();
}



real_t CollisionQuantityHandler::evaluateHColdAtP(len_t i, real_t p) {    
    real_t lnLee;

    // Depending on setting for lnLambda, set to constant or energy dependent function
    if (settings->lnL_type==SimulationGenerator::COLLQTY_LNLAMBDA_CONSTANT)
        lnLee = this->lnLambda_c[i];
    else if (settings->lnL_type==SimulationGenerator::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT)
        lnLee = evaluateLnLambdaEEAtP(i,p);

    // Depending on setting, set nu_s to superthermal or full formula (with maxwellian)
    if (settings->collfreq_mode==SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL)
        return lnLee * constPreFactor * (1+p*p)/(p*p*p);
    else if (settings->collfreq_mode==SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        return -1;// todo: relativistic maxwellian test-particle term
    else
        return -1;
}

real_t CollisionQuantityHandler::evaluateHiAtP(len_t i, real_t p, len_t Z, len_t Z0) {    
    
    if (settings->collfreq_type==SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED)
        if (settings->collfreq_mode==SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL)
            return -1; //todo, bound electrons contribute via Bethe stopping power
        else if (settings->collfreq_mode==SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) 
            return -1; // todo, bound electrons contribute via Bethe stopping power matched to thermal 
                       // nu_s formula with new method for adding nu_|| term (see doc/notes/theory)
        else 
            return -1; // no such setting implemented

    // if not using partial screening models, the definition of n_cold handles whether
    // we have complete or no screening, and the ions therefore do not contribute.
    else 
        return 0;

}

real_t CollisionQuantityHandler::evaluateGColdAtP(len_t i, real_t p) {
    real_t lnLee;

    // Depending on setting for lnLambda, set to constant or energy dependent function
    if (settings->lnL_type==SimulationGenerator::COLLQTY_LNLAMBDA_CONSTANT)
        lnLee = this->lnLambda_c[i];
    else if (settings->lnL_type==SimulationGenerator::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT)
        lnLee = evaluateLnLambdaEEAtP(i,p);

    // Depending on setting, set nu_D to superthermal or full formula (with maxwellian)
    if (settings->collfreq_mode==SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL)
        return lnLee * constPreFactor * sqrt(1+p*p)/(p*p*p);
    else if (settings->collfreq_mode==SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        return -1;// todo: relativistic maxwellian test-particle term
    else
        return -1; // error: no such setting supported

}


real_t CollisionQuantityHandler::evaluateGiAtP(len_t i, real_t p, len_t Z, len_t Z0) {
    real_t lnLei;
    real_t h_i;
    // Depending on setting for lnLambda, set to constant or energy dependent function
    if (settings->lnL_type==SimulationGenerator::COLLQTY_LNLAMBDA_CONSTANT)
        lnLei = this->lnLambda_c[i];
    else if (settings->lnL_type==SimulationGenerator::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT)
        lnLei = evaluateLnLambdaEIAtP(i,p);

    // the completely screened contribution
    h_i = Z0*Z0 * lnLei * constPreFactor * sqrt(1+p*p)/(p*p*p);
    
    if (settings->collfreq_type == SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED)
        h_i += (Z*Z-Z0*Z0) * lnLei * constPreFactor * sqrt(1+p*p)/(p*p*p);
    else if (settings->collfreq_type == SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED)
        h_i = -1; //todo, loads data from file and adds the partially screened contribution
    
    return h_i;
}







real_t CollisionQuantityHandler::evaluateNuSAtP(len_t i, real_t p){
    real_t ns = n_cold[i]*evaluateHColdAtP(i, p);
     for (len_t iZ = 0; iZ<nZ[i]; iZ++)
                    ns += ionDensity[i][iZ]
                        * evaluateHiAtP(i, p, ZAtomicNumber[i][iZ],Z0ChargeNumber[i][iZ]);

    return ns;
                
}
real_t CollisionQuantityHandler::evaluateNuDAtP(len_t i, real_t p){
    real_t nD = n_cold[i]*evaluateGColdAtP(i, p);
     for (len_t iZ = 0; iZ<nZ[i]; iZ++)
                    nD += ionDensity[i][iZ]
                        * evaluateGiAtP(i, p, ZAtomicNumber[i][iZ],Z0ChargeNumber[i][iZ]);
                
    return nD;
}



real_t CollisionQuantityHandler::evaluateLnLambdaC(len_t i) {
    return 14.6 + 0.5*log( T_cold[i]/(n_cold[i]/1e20) );
}

real_t CollisionQuantityHandler::evaluateLnLambdaEEAtP(len_t i, real_t p) {
    real_t gamma = sqrt(p*p+1);
    return evaluateLnLambdaC(i) + log( sqrt(gamma-1) );
}

real_t CollisionQuantityHandler::evaluateLnLambdaEIAtP(len_t i, real_t p) {
    return evaluateLnLambdaC(i) + log( sqrt(2)*p );
}



void CollisionQuantityHandler::CalculateCollisionFrequencies(){
    real_t 
        **nu_s   = new real_t*[n], 
        **nu_D   = new real_t*[n], 
        **nu_D2  = new real_t*[n],
        **nu_s1  = new real_t*[n], 
        **nu_s2  = new real_t*[n], 
        **nu_D1  = new real_t*[n];
    real_t p, p_f1, p_f2;
    bool gridtypePXI, gridtypePPARPPERP;



    for (len_t ir = 0; ir < n; ir++) {
        gridtypePXI         = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP);
        
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        // Terms on distribution grid
        nu_s[ir]  = new real_t[np1*np2];
        nu_D[ir]  = new real_t[np1*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                if (gridtypePXI)
                    p = mg->GetP1(i);
                else if (gridtypePPARPPERP)
                    p = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2(j)*mg->GetP2(j));
                
                nu_s[ir][j*np1+i] = evaluateNuSAtP(ir,p);
                nu_D[ir][j*np1+i] = evaluateNuDAtP(ir,p);
                
            }
        }

        // Terms on flux grid 1
        nu_s1[ir]  = new real_t[(np1+1)*np2];
        nu_D1[ir]  = new real_t[(np1+1)*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                if (gridtypePXI)
                    p_f1 = mg->GetP1_f(i);
                else if (gridtypePPARPPERP)
                    p_f1 = sqrt(mg->GetP1_f(i)*mg->GetP1_f(i) + mg->GetP2(j)*mg->GetP2(j));
                
                nu_s1[ir][j*(np1+1)+i]  = evaluateNuSAtP(ir,p_f1);
                nu_D1[ir][j*(np1+1)+i]  = evaluateNuDAtP(ir,p_f1);
            }
        }

        // Terms on flux grid 2
        nu_s2[ir]  = new real_t[np1*(np2+1)];
        nu_D2[ir]  = new real_t[np1*(np2+1)];
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                if (gridtypePXI)
                    p_f2 = mg->GetP1(i);
                else if (gridtypePPARPPERP)
                    p_f2 = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2_f(j)*mg->GetP2_f(j));
                                
                nu_s2[ir][j*np1+i]  = evaluateNuSAtP(ir,p_f2);
                nu_D2[ir][j*np1+i]  = evaluateNuDAtP(ir,p_f2);
            }
        }
    }
    this->collisionFrequencyNuS    = nu_s;
    this->collisionFrequencyNuS_f1 = nu_s1;
    this->collisionFrequencyNuS_f2 = nu_s2;
    this->collisionFrequencyNuD    = nu_D;
    this->collisionFrequencyNuD_f1 = nu_D1;
    this->collisionFrequencyNuD_f2 = nu_D2;
}


// Loads ADAS coefficients and uses ion species and atomic parmeters data to calculate
//   various ionisation and recombination rates
void CollisionQuantityHandler::CalculateIonisationRates(){

  
    // SetIonisationRates(Icold, Ikin, IRE,
    //                    RR, CEZP, CEHP){
    
}




// Uses collision frequencies and ion species to calculate
// critical fields and avalanche growth rates
void CollisionQuantityHandler::CalculateDerivedQuantities(){

    // SetDerivedQuantities(Ec, Ectot, ED, 
    //                        Gamma, Eceff);
    
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
    }

    delete [] this->collisionFrequencyNuS;
    delete [] this->collisionFrequencyNuS_f1;
    delete [] this->collisionFrequencyNuS_f2;
    delete [] this->collisionFrequencyNuD;
    delete [] this->collisionFrequencyNuD_f1;
    delete [] this->collisionFrequencyNuD_f2;
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

// Bonus function not used by the DREAM simulation workflow
void CollisionQuantityHandler::CalculateCoulombLogarithms(){ 
    if (ionDensity==nullptr){
        // error?
    }
    real_t p, p_f1, p_f2, gamma, gamma_f1, gamma_f2;
    
    real_t 
        **lnLee  = new real_t*[n], 
        **lnLei  = new real_t*[n], 
        **lnLee1 = new real_t*[n], 
        **lnLee2 = new real_t*[n], 
        **lnLei1 = new real_t*[n], 
        **lnLei2 = new real_t*[n],
         *lnLc   = new real_t[n], 
         *lnLTe  = new real_t[n];

    bool gridtypePXI, gridtypePPARPPERP;

    for (len_t ir = 0; ir < n; ir++) {
        gridtypePXI         = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP);
        
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        lnLc[ir]  = 14.6 + 0.5*log( T_cold[ir]/(n_cold[ir]/1e20) );
        lnLTe[ir] = 14.9 + 0.5*log( (T_cold[ir]/1e3)*(T_cold[ir]/1e3)/(n_cold[ir]/1e20) );

        lnLee[ir] = new real_t[np1*np2];
        lnLei[ir] = new real_t[np1*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                if (gridtypePXI)
                    p = mg->GetP1(i);
                else if (gridtypePPARPPERP)
                    p = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2(j)*mg->GetP2(j));
                gamma = sqrt(1+p*p);
                lnLee[ir][j*np1+i] = lnLc[ir] + log( sqrt(gamma-1) );
                lnLei[ir][j*np1+i] = lnLc[ir] + log( sqrt(2)*p );
            }
        }

        lnLee1[ir] = new real_t[(np1+1)*np2];
        lnLei1[ir] = new real_t[(np1+1)*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                if (gridtypePXI)
                    p_f1 = mg->GetP1_f(i);
                else if (gridtypePPARPPERP)
                    p_f1 = sqrt(mg->GetP1_f(i)*mg->GetP1_f(i) + mg->GetP2(j)*mg->GetP2(j));
                gamma_f1 = sqrt(1+p_f1*p_f1);
                lnLee1[ir][j*(np1+1)+i] = lnLc[ir] + log( sqrt(gamma_f1-1) );
                lnLei1[ir][j*(np1+1)+i] = lnLc[ir] + log( sqrt(2)*p_f1 );
            }
        }

        lnLee2[ir] = new real_t[np1*(np2+1)];
        lnLei2[ir] = new real_t[np1*(np2+1)];
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                if (gridtypePXI)
                    p_f2 = mg->GetP1(i);
                else if (gridtypePPARPPERP)
                    p_f2 = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2_f(j)*mg->GetP2_f(j));
                gamma_f2 = sqrt(1+p_f2*p_f2);
                lnLee2[ir][j*np1+i] = lnLc[ir] + log( sqrt(gamma_f2-1) );
                lnLei2[ir][j*np1+i] = lnLc[ir] + log( sqrt(2)*p_f2 );
            }
        }
    }


}







