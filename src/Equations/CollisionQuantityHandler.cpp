/**
 * So far this is the outline of a class that could calculate
 * various collision, ionisation, ... rates and quantities. 
 * You would set ion species via a 
 * list of ion densities, described by atomic numbers Z and 
 * charge numbers Z0. Then it loads atomic physics data and
 * calculates a bunch of collison-related quantities.
 * Could also evaluate avalanche growth rates and critical fields.
*/


/** 
 * EXAMPLE:
 * 
 * /---/
 * CollisionQuantityHandler *CollQty;

 * CollQty->SetIonSpecies(densities, Zs, Z0s)
 * 
 * real_t *lnLee = CollQty->evaluateLnLeeAtP(p);
 * real_t *nu_s  = CollQty->evaluateNuSAtP(p);
 * 
 * CollQty->SetGrid(grid);
 * CollQty->CalculateCollisionFrequencies(); // stores nu_s and nu_D on full grid
 * real_t **nu_s = CollQty->GetNuS();
 * 
 * CollQty->SetSpeciesFromEqSys(equationSystem);
 * CollQty->CalculateIonisationRates();
 */



#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

#include <cmath>





using namespace DREAM;


/** 
 * Constructor
 */ 
CollisionQuantityHandler::CollisionQuantityHandler(struct collqtyhand_settings *cq){
    // ?
}


/**
 * Destructor.
 */
CollisionQuantityHandler::~CollisionQuantityHandler(){
    DeallocateCollisionFrequencies();
    DeallocateIonisationRates();
    DeallocateIonSpecies();
}


void CollisionQuantityHandler::DeallocateCollisionFrequencies(){
    if (this->lnLambda_c == nullptr)
        return;

    delete [] this->lnLambda_c;
    delete [] this->lnLambda_Te;
    

    for (len_t i = 0; i < this->n; i++) {
        delete [] this->lnLambda_ee[i];
        delete [] this->lnLambda_ei[i];
        delete [] this->collisionFrequencyNuS[i];
        delete [] this->collisionFrequencyNuD[i];
    }

    delete [] this->lnLambda_ee;
    delete [] this->lnLambda_ei;
    delete [] this->collisionFrequencyNuS;
    delete [] this->collisionFrequencyNuD;
}

/**
 * Uses ion species and atomic parameters data to calculate
 * collision frequencies. How to choose settings?
 */

void CollisionQuantityHandler::CalculateCollisionFrequencies(){

    /*
    ...

     SetCollisionFrequencies(nu_s, nu_D, lnLc,
                                lnLee, lnLei, lnLTe);
    */
}


// Loads ADAS coefficients and uses ion species and atomic parmeters data to calculate
// various ionisation and recombination rates
void CollisionQuantityHandler::CalculateIonisationRates(){

    /*
    ...

    SetIonisationRates(Icold, Ikin, IRE,
                        RR, CEZP, CEHP){
    */
}

// Uses collision frequencies and ion species to calculate
// critical fields and avalanche growth rates
void CollisionQuantityHandler::CalculateDerivedQuantities(){

    /*
    ...

    SetDerivedQuantities(Ec, Ectot, ED, 
                            Gamma, Eceff);
    */
}











