/**
 * Implementation of collision-rate calculator that calculates
 * various collision, ionisation etc rates and quantities.  
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
#include "DREAM/EquationSystem.hpp"

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



void CollisionQuantityHandler::RebuildFromEqSys() {
    /** 
     * 
     * *extract quantities from EqSys* (all densities and charge states)
     * 
     * SetIonSpecies(ionDensities, Zs, Z0s);
     * 
     * this->T_cold = T;
     * this->n_cold = n_cold;
     * this->n_fast = n_hot + n_RE;
     * 
     * CalculateLnLambda();        // since temperatures have changed
     * CalculateIonisationRates(); // like above
     * 
     * if (grid or Zs or Z0s changed)
     *     LoadAtomicData()
     *     CalculateGiFunctions(); // calculates functions g_i... when to deallocate? pretty significant memory cost
     *                             // as they live on nZ x (n1 x n2) grid. Of course, could optimize for
     *                             // PXI grid where they could live on nZ x n1.
     *                             // Could probably deallocate once jacobian is built. 
     * CalculateCollisionFrequencies(); // nu_s = n_cold*lnLambda_ee*g_{s,e} + sum_i n_i g_{s,i}
     *                                  // nu_D = n_cold*lnLambda_ei*g_{D,e} + sum_i n_i g_{D,i}
     * CalculateDerivedQuantities(); 
     * 
     */


}

// Revise the below: 
void CollisionQuantityHandler::SetIonSpecies(real_t **dens, len_t **Z, len_t **Z0){
    DeallocateIonSpecies();
    this->ionDensity     = dens;
    this->ZAtomicNumber  = Z;
    this->Z0ChargeNumber = Z0;



    if (eqSys == nullptr) {
        real_t *nZeff  = new real_t[n];
        real_t *nZeff0 = new real_t[n];
        real_t *nfree  = new real_t[n];
        real_t *ntot   = new real_t[n];
        
        for (len_t i   = 0; i<n; i++){
            nfree[i]   = 0;
            ntot[i]    = 0;
            nZeff[i]   = 0;
            nZeff0[i]  = 0;
            for (len_t iz = 0; iz<nZ; iz++){
                nfree[i]  += Z0[i][iz]*dens[i][iz];
                ntot[i]   += Z[i][iz]*dens[i][iz];
                nZeff0[i] += Z0[i][iz]*Z0[i][iz]*dens[i][iz];
                nZeff[i]  += Z[i][iz]*Z[i][iz]*dens[i][iz];
            }
        }

        this->n_free = nfree;
        this->n_total = ntot;
        this->nTimesZeff0ForNuD = nZeff0;
        this->nTimesZeffForNuD = nZeff;
    }
    
    /** 
     * Right now: if selfconsistent it is assumed that SetSpeciesFromEqSys 
     * has been run which sets n_fast and n_cold.
     *
     * Actually, should the simulation really do anything with n_cold? I think
     * that the only place it is actually used is in this class. It is not a quantity 
     * that we solve for, but just obtained via impurities, n_hot and n_RE via 
     * quasi charge neutrality.
     *
     *
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
    */
    
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
        *lnLc    = new real_t[n], 
        *lnLTe   = new real_t[n];

    bool gridtypePXI, gridtypePPARPPERP;

    // first add contribution from free electrons
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




/**
 * CHANGE: the below will be revised by splitting the momentum-space 
 * functions and the densities, so that the calculation will look like
 * nu_s = nu_sfree + \sum_i n_i g_i
 * nu_D = nu_sfree + \sum_i n_i g_i
 * where the sum is taken over all ion species in the plasma
 * and the functions g_i are calculated in a separate function (and stored? at least until 
 * the Jacobian matrix for the newton solver has been assembled).
 * Regardless, it will probably be a good idea to define a function CalculateGi(ir, p, Z, Z0)
 * which would take care of much of the code repeated below.
 */
/**
 * Uses ion species and atomic parameters data to calculate
 * collision frequencies. 
 */
void CollisionQuantityHandler::CalculateCollisionFrequencies(){
    if (ionDensity==nullptr){
        // error?
    }
    real_t p, p_f1, p_f2, gamma, gamma_f1, gamma_f2;
    

    real_t *n_nu_s, *Zeff_n_nu_D; // the density of free electrons that we use in nu_s and nu_D prefactors
    if (settings->collfreq_type == SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED) {
        n_nu_s      = this->n_total;
        Zeff_n_nu_D = this->nTimesZeffForNuD;
    }
    else if (settings->collfreq_type == SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED) {
        n_nu_s      = this->n_cold;
        Zeff_n_nu_D = this->nTimesZeff0ForNuD;
    }


    real_t 
        **nu_s   = new real_t*[n], 
        **nu_D   = new real_t*[n], 
        **nu_D2  = new real_t*[n],
        **nu_s1  = new real_t*[n], 
        **nu_s2  = new real_t*[n], 
        **nu_D1  = new real_t*[n],
        **commonMomentumFunction    = new real_t*[n],
        **commonMomentumFunction_f1 = new real_t*[n],
        **commonMomentumFunction_f2 = new real_t*[n];

    real_t  
        **lnLee  = GetLnLambdaEE(),
        **lnLei  = GetLnLambdaEI(),
        **lnLee1 = GetLnLambdaEE_f1(), 
        **lnLee2 = GetLnLambdaEE_f2(),
        **lnLei1 = GetLnLambdaEI_f1(),
        **lnLei2 = GetLnLambdaEI_f2(),
         *lnLc   = GetLnLambdaC(),
         *lnLTe  = GetLnLambdaTe();
        
        
    bool gridtypePXI, gridtypePPARPPERP;

    // first add contribution from free electrons
    for (len_t ir = 0; ir < n; ir++) {
        gridtypePXI         = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP);
        
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        nu_s[ir]  = new real_t[np1*np2];
        nu_D[ir]  = new real_t[np1*np2];
        commonMomentumFunction[ir]    = new real_t[np1*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                if (gridtypePXI)
                    p = mg->GetP1(i);
                else if (gridtypePPARPPERP)
                    p = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2(j)*mg->GetP2(j));
                gamma = sqrt(1+p*p);
                
                commonMomentumFunction[ir][j*np1+i] = constPreFactor*gamma/(p*p*p);
                
                nu_s[ir][j*np1+i]  = n_nu_s[ir] * lnLee[ir][j*np1+i]
                                    * gamma*commonMomentumFunction[ir][j*np1+i];
                nu_D[ir][j*np1+i]  =  Zeff_n_nu_D[ir] * constPreFactor 
                                    * lnLee[ir][j*np1+i] * commonMomentumFunction[ir][j*np1+i];
            }
        }

        nu_s1[ir]  = new real_t[(np1+1)*np2];
        nu_D1[ir]  = new real_t[(np1+1)*np2];
        commonMomentumFunction_f1[ir] = new real_t[(np1+1)*np2];
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                if (gridtypePXI)
                    p = mg->GetP1_f(i);
                else if (gridtypePPARPPERP)
                    p = sqrt(mg->GetP1_f(i)*mg->GetP1_f(i) + mg->GetP2(j)*mg->GetP2(j));
                gamma = sqrt(1+p*p);
                
                commonMomentumFunction_f1[ir][j*(np1+1)+i] = constPreFactor*gamma/(p*p*p);
                
                nu_s1[ir][j*(np1+1)+i]  = n_nu_s[ir] * lnLee1[ir][j*(np1+1)+i]
                                    * gamma*commonMomentumFunction_f1[ir][j*(np1+1)+i];
                nu_D1[ir][j*(np1+1)+i]  = Zeff_n_nu_D[ir] * constPreFactor 
                                    * lnLee1[ir][j*(np1+1)+i] * commonMomentumFunction_f1[ir][j*(np1+1)+i];
            }
        }

        nu_s2[ir]  = new real_t[np1*(np2+1)];
        nu_D2[ir]  = new real_t[np1*(np2+1)];
        commonMomentumFunction_f2[ir]    = new real_t[np1*(np2+1)];
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                if (gridtypePXI)
                    p = mg->GetP1(i);
                else if (gridtypePPARPPERP)
                    p = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2_f(j)*mg->GetP2_f(j));
                gamma = sqrt(1+p*p);
                
                commonMomentumFunction_f2[ir][j*np1+i] = constPreFactor*gamma/(p*p*p);
                
                nu_s2[ir][j*np1+i]  = n_nu_s[ir] * lnLee2[ir][j*np1+i]
                                    * gamma*commonMomentumFunction_f2[ir][j*np1+i];
                nu_D2[ir][j*np1+i]  = Zeff_n_nu_D[ir] * constPreFactor 
                                    * lnLee2[ir][j*np1+i] * commonMomentumFunction_f2[ir][j*np1+i];
            }
        }
    }



    // If using partial screening model, sum over ion species and weigh with atomic data
    if (settings->collfreq_type == SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED) {
        
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










