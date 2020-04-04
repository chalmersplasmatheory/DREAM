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
            enum SimulationGenerator::collqty_collfreq_type collfreq_type   = SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED;
            enum SimulationGenerator::collqty_collfreq_mode collfreq_mode   = SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
            enum SimulationGenerator::collqty_lnLambda_type lnL_type        = SimulationGenerator::COLLQTY_LNLAMBDA_CONSTANT;
            enum SimulationGenerator::uqty_n_cold_eqn       ncold_type      = SimulationGenerator::UQTY_N_COLD_EQN_PRESCRIBED;
        };

    private:
        const real_t constPreFactor = 4*Constants::pi 
                                *Constants::r0*Constants::r0
                                *Constants::c;
        len_t n,   // number of "radial grid points" (or sets of ion species) 
              *nZ; // number of ion species at n
        FVM::Grid *grid;
        EquationSystem *eqSys = nullptr;
        enum SimulationGenerator::momentumgrid_type gridtype;

        // Ion densities on n x nZ
        real_t  *n_cold;                 // thermal free electron density in m^-3
        real_t  *T_cold;                 // thermal free electron temperature in eV
        real_t **ionDensity=nullptr;     // ion densities in m^-3
        len_t  **ZAtomicNumber;          // atomic number (nuclear charge) of ion
        len_t  **Z0ChargeNumber;         // charge number (net charge) of ion

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
        real_t **collisionFrequencyNuD_f2=nullptr; // pitch scatter frequency on flux grid 2
        

        // Partially-screened ion contribution to nu_s and nu_D on n x nZ x (np1 x np2)
        real_t ***HiPartialScreened;
        real_t ***HiPartialScreened_f1;
        real_t ***HiPartialScreened_f2;
        real_t ***GiPartialScreened;
        real_t ***GiPartialScreened_f1;
        real_t ***GiPartialScreened_f2=nullptr;

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


        void SetEqSys(EquationSystem *es){
            this->eqSys = es;
        }
        virtual void RebuildFromEqSys();

        /** 
         * The g_i and h_i functions are defined so that 
         * nu_s = n_cold*h_cold + sum_i n_i*h_i, 
         * nu_D = n_cold*g_cold + sum_i n_i*g_i,
         * where i sums over all different ion species (ie combinations of Z and Z0)
         * and g_i is purely contribution from _bound_ electrons
         * 
         * Currently we don't store the partially-screened contributions to g_i and h_i, which we probably
         * should since they are always the same (but somewhat costly with memory since it is on an 
         * nZ x np1 x np2 array) 
         */
        virtual real_t evaluateGiAtP(len_t i, real_t p, len_t Z, len_t Z0);
        virtual real_t evaluateHiAtP(len_t i, real_t p, len_t Z, len_t Z0);
        virtual real_t evaluateGColdAtP(len_t i, real_t p);
        virtual real_t evaluateHColdAtP(len_t i, real_t p);
        
        virtual real_t evaluateKirillovGiAtP(real_t p, len_t Z, len_t Z0);
        virtual real_t evaluateBetheHiAtP(real_t p, len_t Z, len_t Z0);
        
        virtual real_t evaluateNuSAtP(len_t i, real_t p);
        virtual real_t evaluateNuDAtP(len_t i, real_t p);
        
        real_t evaluateLnLambdaEEAtP(len_t i,real_t p);
        real_t evaluateLnLambdaEIAtP(len_t i,real_t p);
        real_t evaluateLnLambdaC(len_t i);

        virtual void CalculateCollisionFrequencies(); // lnL and nu

        virtual void CalculateIonisationRates();      // I, R and CE
        virtual void CalculateDerivedQuantities();    // Ec, Gamma_ava
        


        virtual void DeallocateCollisionFrequencies();
        virtual void DeallocateIonisationRates();
        virtual void DeallocateIonSpecies();
        virtual void DeallocateLnLambdas();
        virtual void DeallocateDerivedQuantities();
        virtual void DeallocateHiGiPartialScreened();
        virtual real_t GetIonEffectiveSizeAj(len_t Z, len_t Z0);   // search atomic-data table for the Z, Z0 value. 
        virtual real_t GetMeanExcitationEnergy(len_t Z, len_t Z0); // search atomic-data table for the Z, Z0 value



        real_t *const* GetNuS() const 
                { return this->collisionFrequencyNuS; }
        const real_t  *GetNuS(const len_t i) const 
                { return this->collisionFrequencyNuS[i]; }
        real_t *const* GetNuS_f1() const 
                { return this->collisionFrequencyNuS_f1; }
        const real_t  *GetNuS_f1(const len_t i) const 
                { return this->collisionFrequencyNuS_f1[i]; }
        real_t *const* GetNuS_f2() const 
                { return this->collisionFrequencyNuS_f2; }
        const real_t  *GetNuS_f2(const len_t i) const 
                { return this->collisionFrequencyNuS_f2[i]; }

        real_t *const* GetNuD() const 
                { return this->collisionFrequencyNuD; }
        const real_t  *GetNuD(const len_t i) const 
                { return this->collisionFrequencyNuD[i]; }
        real_t *const* GetNuD_f1() const 
                { return this->collisionFrequencyNuD_f1; }
        const real_t  *GetNuD_f1(const len_t i) const 
                { return this->collisionFrequencyNuD_f1[i]; }
        real_t *const* GetNuD_f2() const 
                { return this->collisionFrequencyNuD_f2; }
        const real_t  *GetNuD_f2(const len_t i) const 
                { return this->collisionFrequencyNuD_f2[i]; }




        // this one may be used if we want to sacrifice memory for cpu performance?
        virtual void CalculateHiGiPartialScreened(); 
        
        /**
         * NOTE: The below methods are not used in the standard DREAM workflow
         */


        /**
         * Calculates and stores g and h functions on grid,
         * on nZ x n x (np1 x np2) arrays. 
         * The benefit of implementing this method is that we need to 
         * calculate all these functions twice: once for the advection
         * and diffusion terms, and once for the Jacobian matrix for 
         * the newton solver.
         */ 
        virtual void CalculateGiHiFuncs();
        
        virtual void CalculateCoulombLogarithms(); // lnL and nu
                
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
        
        real_t *const* GetIonDens() const 
                { return this->ionDensity; }
        const real_t  *GetIonDens(const len_t i) const 
                { return this->ionDensity[i]; }
        const real_t  GetIonDens(const len_t i, const len_t Z, const len_t Z0) const { 
            for(len_t iz = 0; iz<nZ[i]; iz++) {
                if ( ZAtomicNumber[i][iz]==Z && Z0ChargeNumber[i][iz]==Z0 )
                    return ionDensity[i][iz];
            }
        }
        // and so on 


        void SetGrid(FVM::Grid *g, enum SimulationGenerator::momentumgrid_type mgtype){
            this->grid = g;
            this->gridtype = mgtype;

        }


        void SetTemperature(real_t *T){
            this->T_cold = T;
        }


        virtual void SetIonSpecies(real_t **dens, len_t **Z, len_t **Z0, real_t *T);

        // is this needed?
        void SetCollisionFrequencies(
                real_t **nu_s, real_t **nu_D, real_t **nu_D2,
                real_t **nu_s1, real_t **nu_s2, real_t **nu_D1, 
                real_t *lnLc, real_t **lnLee, real_t **lnLei, 
                real_t *lnLTe, real_t **lnLee1, real_t **lnLee2, 
                real_t **lnLei1, real_t **lnLei2
            ) {
            DeallocateCollisionFrequencies();
            DeallocateLnLambdas();
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


        
        virtual void LoadAtomicData();

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


      
        
        
    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP*/

    


