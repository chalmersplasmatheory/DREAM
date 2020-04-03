#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP


#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM {
    class CollisionQuantityHandler{

    public:
        struct collqtyhand_settings {
            enum SimulationGenerator::collqty_collfreq_type collfreq_type=SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_SUPERTHERMAL_NON_SCREENED;
            enum SimulationGenerator::uqty_n_cold_eqn ncold_type = SimulationGenerator::UQTY_N_COLD_EQN_PRESCRIBED;
        };

    private:
        const real_t constPreFactor = 4*Constants::pi 
                                *Constants::r0*Constants::r0
                                *Constants::c;
        len_t n,  // number of "radial grid points" (or sets of ion species) 
              nZ; // number of ion species
        FVM::Grid *grid;
        EquationSystem *eqSys = nullptr;
        enum SimulationGenerator::momentumgrid_type gridtype;
        
        
        // For each "radial position" (or just index) i < n
        // you can have a number nZ of different ion species 
        // (i.e. deuterium+1, neon+3, argon+0, argon+1 => nZ=4)  

        // Ion densities on n x nZ
        real_t  *n_cold;                 // thermal free electron density in m^-3
        real_t  *n_free;                 // total free electron density in m^-3
        real_t  *n_total;                // total free+bound electron density in m^-3
        real_t  *nTimesZeff0ForNuD;       // n*Zeff as appearing in completely screened nuD
        real_t  *nTimesZeffForNuD;       // n*Zeff as appearing in non-screened nuD
        real_t  *n_fast;                 // equals n_hot + n_RE in m^-3
        real_t  *T_cold;                 // thermal free electron temperature in eV
        real_t **ionDensity=nullptr;     // ion densities in m^-3
        len_t  **ZAtomicNumber;          // atomic number (nuclear charge) of ion
        len_t  **Z0ChargeNumber;         // charge number (net charge) of ion
        len_t   *len_Zvec;
        real_t  *Zeff;                   // effective charge defined such that the 
                                         // completely screened nu_D is propto ncold*Zeff;
        // Coulomb logarithms on n x (np1 x np2)
        real_t  *lnLambda_c=nullptr;     // constant relativistic lnLambda
        real_t **lnLambda_ee;            // energy dependent ee lnLambda
        real_t **lnLambda_ee_f1;         // energy dependent ee lnLambda on flux grid 1
        real_t **lnLambda_ee_f2;         // energy dependent ee lnLambda on flux grid 2
        real_t **lnLambda_ei;            // energy dependent ei lnLambda
        real_t **lnLambda_ei_f1;         // energy dependent ei lnLambda on flux grid 1
        real_t **lnLambda_ei_f2;         // energy dependent ei lnLambda on flux grid 2
        real_t  *lnLambda_Te;            // constant thermal lnLambda
        

        // Collision frequencies on n x (np1 x np2)   
        real_t **collisionFrequencyNuS;    // slowing down frequency
        real_t **collisionFrequencyNuS_f1; // slowing down frequency on flux grid 1
        real_t **collisionFrequencyNuS_f2; // slowing down frequency on flux grid 2
        real_t **collisionFrequencyNuD;    // pitch scatter frequency
        real_t **collisionFrequencyNuD_f1; // pitch scatter frequency on flux grid 1
        real_t **collisionFrequencyNuD_f2; // pitch scatter frequency on flux grid 2
        
        // Ionisation and recombination on n x nZ
        real_t **ionisationRateCold=nullptr;  // ionisation rate by thermal cold electrons
        real_t **ionisationRateKinetic;       // ionisation rate by kinetic electrons
        //         ^-- perhaps this should be a function of p? (or even f?)
        //             unreasonable to store on n x nZ x (np1 x np2)..?
        real_t **ionisationRateRE;            // ionisation rate by RE fluid
        real_t **recombinationRateRadiative;  // radiative recombination rate
        real_t **chargeExchangeZP;            // impurity-proton charge exchange rates
        real_t  *chargeExchangeHP;            // hydrogen-proton charge exchange rate
        
        
        // Atomic parameters on nZ
        real_t *ionisationPotential=nullptr;  // Ionisation energies loaded from file
        real_t *meanExcitationEnergy;         // For nu_s. Loaded from file or calculated.
        real_t *ionEffectiveSizeAj;           // For nu_D. Loaded from file or calculated.

        // Kinetic derived quantities on n grid
        real_t *Ec_free=nullptr;        // Connor-Hastie field with only bound
        real_t *Ec_tot;                 // Connor-Hastie field with free+bound
        real_t *EDreic;                 // Dreicer field
        real_t *avalancheGrowthRate;    // Gamma_ava
        real_t *effectiveCriticalField; // Eceff: Gamma_ava(Eceff) = 0

        struct collqtyhand_settings *settings;

    public:

        CollisionQuantityHandler(struct collqtyhand_settings *cq=nullptr);
        ~CollisionQuantityHandler();

        // once ion species have been set, various 
        // quantities can be calculated. 
        virtual void CalculateCollisionFrequencies(); // lnL and nu
        virtual void CalculateCoulombLogarithms(); // lnL and nu
        
        /**
         * lnL and nu matched to screened thermal collision rates for p≲p_Te.
         * Can probably also create nu_|| in such a way as to make
         * the relativistic maxwellian f_MR an exact solution, i.e. by writing the collisional p-flux
         * as 1/p^2 d/dp [ p^4*nu_col* f_MR d/dp(f/f_MR) ]   
         * where nu_col should be chosen so as to make the total energy loss rate equal to the stopping power ("vpnu_s")
         */
        virtual void CalculateCollisionFrequenciesWithFiniteT(); 
        
        virtual void CalculateIonisationRates();      // I, R and CE
        virtual void CalculateDerivedQuantities();    // Ec, Gamma_ava
         
        real_t *evaluateNuSAtP(real_t p);
        // and so on

        
        real_t **GetLnLambdaEE()  
                { return this->lnLambda_ee; }
        real_t **GetLnLambdaEI()  
                { return this->lnLambda_ei; }
        real_t **GetLnLambdaEE_f1()  
                { return this->lnLambda_ee_f1; }
        real_t **GetLnLambdaEE_f2()  
                { return this->lnLambda_ee_f2; }
        real_t **GetLnLambdaEI_f1()  
                { return this->lnLambda_ei_f1; }
        real_t **GetLnLambdaEI_f2()  
                { return this->lnLambda_ei_f2; }
        real_t *GetLnLambdaC() 
                { return this->lnLambda_c; }
        real_t *GetLnLambdaTe() 
                { return this->lnLambda_Te; }
        
        /* kunde inte använda detta som jag önskade i CalculateCollisionFrequencies..?
        real_t *const* GetLnLambdaEE() const 
                { return this->lnLambda_ee; }
        real_t *const* GetLnLambdaEI() const 
                { return this->lnLambda_ei; }
        real_t *const* GetLnLambdaEE_f1() const 
                { return this->lnLambda_ee_f1; }
        real_t *const* GetLnLambdaEE_f2() const 
                { return this->lnLambda_ee_f2; }
        real_t *const* GetLnLambdaEI_f1() const 
                { return this->lnLambda_ei_f1; }
        real_t *const* GetLnLambdaEI_f2() const 
                { return this->lnLambda_ei_f2; }
        const real_t *GetLnLambdaC() const
                { return this->lnLambda_c; }
        const real_t *GetLnLambdaTe() const
                { return this->lnLambda_Te; }
        */
        real_t *const* GetNuS() const 
                { return this->collisionFrequencyNuS; }
        const real_t  *GetNuS(const len_t i) const 
                { return this->collisionFrequencyNuS[i]; }
        real_t *const* GetIonDens() const 
                { return this->ionDensity; }
        const real_t  *GetIonDens(const len_t i) const 
                { return this->ionDensity[i]; }
        const real_t  GetIonDens(const len_t i, const len_t Z, const len_t Z0) const { 
            for(len_t n = 0; n<len_Zvec[i]; n++) {
                if ( ZAtomicNumber[i][n]==Z && Z0ChargeNumber[i][n]==Z0 )
                    return ionDensity[i][n];
            }
        }
        // and so on 


        void SetGrid(FVM::Grid *g, enum SimulationGenerator::momentumgrid_type mgtype){
            this->grid = g;
            this->gridtype = mgtype;

        }

        void SetEqSys(EquationSystem *es){
            this->eqSys = es;
        }

        virtual void RebuildFromEqSys();

        virtual void SetIonSpecies(real_t **dens, len_t **Z, len_t **Z0);

        // is the idea that the user should be able to set custom models in a simulation via this method?
        // in that case revise with nu_s = sum_i n_i g_i treatment.
        void SetCollisionFrequencies(
                real_t **nu_s, real_t **nu_D, real_t **nu_D2,
                real_t **nu_s1, real_t **nu_s2, real_t **nu_D1, 
                real_t *lnLc, real_t **lnLee, real_t **lnLei, 
                real_t *lnLTe, real_t **lnLee1, real_t **lnLee2, 
                real_t **lnLei1, real_t **lnLei2
            ) {
            DeallocateCollisionFrequencies();
            this->collisionFrequencyNuS    = nu_s;
            this->collisionFrequencyNuS_f1 = nu_s1;
            this->collisionFrequencyNuS_f2 = nu_s2;
            this->collisionFrequencyNuD    = nu_D;
            this->collisionFrequencyNuD_f1 = nu_D1;
            this->collisionFrequencyNuD_f2 = nu_D2;
            this->lnLambda_c     = lnLc;
            this->lnLambda_ee    = lnLee;
            this->lnLambda_ee_f1 = lnLee1;
            this->lnLambda_ee_f2 = lnLee2;
            this->lnLambda_ei    = lnLei;
            this->lnLambda_ei_f1 = lnLei1;
            this->lnLambda_ei_f2 = lnLei2;
            this->lnLambda_Te    = lnLTe;
        }

        void SetAtomicParameters(real_t *I, real_t *Imean){
            this->meanExcitationEnergy = Imean;
            this->ionisationPotential = I;
        }


        void SetIonisationRates(real_t **Icold, real_t **Ikin, real_t **IRE,
                                    real_t **RR, real_t **CEZP, real_t *CEHP){
            DeallocateIonisationRates();
            this->ionisationRateCold = Icold;
            this->ionisationRateKinetic = Ikin;
            this->ionisationRateRE = IRE;
            this->recombinationRateRadiative = RR;
            this->chargeExchangeZP = CEZP;
            this->chargeExchangeHP = CEHP;       
        }


        void SetDerivedQuantities(real_t *Ec, real_t *Ectot, real_t *ED, 
                                                    real_t *Gamma,  real_t *Eceff){ 
            this->Ec_free = Ec;
            this->Ec_tot  =  Ectot;
            this->EDreic  =  ED;
            this->avalancheGrowthRate = Gamma;
            this->effectiveCriticalField = Eceff;
        }


      
        virtual void DeallocateCollisionFrequencies();
        virtual void DeallocateIonisationRates();
        virtual void DeallocateIonSpecies();
        
    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP*/

    


