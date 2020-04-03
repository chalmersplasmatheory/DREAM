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
    DeallocateIonSpecies();
}

void CollisionQuantityHandler::SetIonSpecies(real_t **dens, len_t **Z, len_t **Z0){
    DeallocateIonSpecies();
    this->ionDensity     = dens;
    this->ZAtomicNumber  = Z;
    this->Z0ChargeNumber = Z0;


    /**
     * this loop defined n_cold (and Zeff) if n_cold is prescribed 
     * (ie determined from ion densities). Otherwise nc below is 
     * n_cold + n_hot + n_RE
     */
    
    real_t *Ze     = new real_t[n];
    real_t *nfree = new real_t[n];  
    for (len_t i = 0; i<n; i++){
        nfree[i] = 0;
        Ze[i] = 0;
        for (len_t iz = 0; iz<nZ; iz++){
            nfree[i] += Z0[i][iz]*dens[i][iz];
            Ze[i]    += Z0[i][iz]*Z0[i][iz]*dens[i][iz];
        }
    }
    
    /** 
     * Right now: if selfconsistent it is assumed that SetSpeciesFromEqSys 
     * has been run which sets n_fast and n_cold
     */
    
    if (settings->ncold_type == SimulationGenerator::UQTY_N_COLD_EQN_PRESCRIBED) {
        for (len_t i = 0; i<n; i++)
            Ze[i] = Ze[i] / nfree[i];
                
        this->Zeff   = Ze;
        this->n_cold = nfree;
        this->n_fast = 0;

    } else if (settings->ncold_type == SimulationGenerator::UQTY_N_COLD_EQN_SELFCONSISTENT) {
        for (len_t i = 0; i<n; i++)
            Ze[i] = Ze[i] / this->n_cold[i];

        this->Zeff = Ze;  
    }
    
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
    if (ionDensity==nullptr){
        // error?
    }
    real_t p, p_f1, p_f2, gamma, gamma_f1, gamma_f2;
    // first add contribution from free electrons
    
    real_t 
        **nu_s   = new real_t*[n], 
        **nu_D   = new real_t*[n], 
        **nu_D2  = new real_t*[n],
        **nu_s1  = new real_t*[n], 
        **nu_s2  = new real_t*[n], 
        **nu_D1  = new real_t*[n], 
        **lnLee  = new real_t*[n], 
        **lnLei  = new real_t*[n], 
        **lnLee1 = new real_t*[n], 
        **lnLee2 = new real_t*[n], 
        **lnLei1 = new real_t*[n], 
        **lnLei2 = new real_t*[n],
        *lnLc    = new real_t[n], 
        *lnLTe   = new real_t[n];

    bool gridtypePXI, gridtypePPARPPERP;                         
    for (len_t ir = 0; ir < n; ir++) {
        gridtypePXI         = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP);
        
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        lnLc[ir]  = 14.6 + 0.5*log( T_cold[ir]/(n_cold[ir]/1e20) );
        lnLTe[ir] = 14.9 + 0.5*log( (T_cold[ir]/1e3)*(T_cold[ir]/1e3)/(n_cold[ir]/1e20) );


        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                if (gridtypePXI)
                    p = mg->GetP1(i);
                else if (gridtypePPARPPERP)
                    p = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2(j)*mg->GetP2(j));
                gamma = sqrt(1+p*p);

                lnLee[ir][j*np1+i] = lnLc[ir] + log( sqrt(gamma-1) );
                lnLei[ir][j*np1+i] = lnLc[ir] + log( sqrt(2)*p );
                nu_s[ir][j*np1+i]  = constPreFactor * n_cold[ir] * lnLee[ir][j*np1+i]
                                    * gamma*gamma/(p*p*p);
                nu_D[ir][j*np1+i]  = constPreFactor * (1+Zeff[ir]) * n_cold[ir]
                                    * lnLee[ir][j*np1+i] * gamma/(p*p*p);

            }
        }

    }

    switch (settings->collfreq_type) {
        case SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED:
            // do nothing, probably?
            break;
        case SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED:
            // add contribution from bound electrons as if they were free 
            break;
        case SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED:
            // sum over not-fully-ionised ion densities with mean-excitation-energy factor
            break;
        default:
            //error: invalid collqty_collfreq_type ?
            ;
    }
 

    /*
    ...

     SetCollisionFrequencies(nu_s, nu_D, nu_D2, nu_s1, nu_s2, nu_D1, 
                lnLc, lnLee, lnLei, lnLTe, lnLee1, lnLee2, lnLei1, lnLei2);
    */
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










