#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP


#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Constants.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>

namespace DREAM {
    class CollisionQuantityHandler{

    public:
        struct collqtyhand_settings {
            enum SimulationGenerator::collqty_collfreq_type 
                        collfreq_type   = SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED;
            enum SimulationGenerator::collqty_collfreq_mode 
                        collfreq_mode   = SimulationGenerator::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
            enum SimulationGenerator::collqty_lnLambda_type 
                        lnL_type        = SimulationGenerator::COLLQTY_LNLAMBDA_CONSTANT;
            enum SimulationGenerator::uqty_n_cold_eqn       
                        ncold_type      = SimulationGenerator::UQTY_N_COLD_EQN_PRESCRIBED;
        };

    private:
        const real_t constPreFactor = 4*M_PI
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
        
        real_t **collisionFrequencyNuPar;
        real_t **collisionFrequencyNuPar_f1;
        real_t **collisionFrequencyNuPar_f2;

        // Partially-screened ion contribution to nu_s and nu_D 
        // on n x (np1 x np2) x nZ
        real_t ***HiFunc;
        real_t ***HiFunc_f1;
        real_t ***HiFunc_f2;
        real_t ***GiFunc;
        real_t ***GiFunc_f1;
        real_t ***GiFunc_f2=nullptr;
        real_t **HCold;
        real_t **HCold_f1;
        real_t **HCold_f2;
        real_t **GCold;
        real_t **GCold_f1;
        real_t **GCold_f2;

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

        
        // atomic data in no particular order, but ...data[i] corresponds to the value for charge Z = ...Zs[i] and Z0 = ...Z0s[i]
        static const len_t ionSizeAj_len; 
        static const real_t ionSizeAj_data[];       
        static const len_t ionSizeAj_Zs[];
        static const len_t ionSizeAj_Z0s[];

        static const len_t meanExcI_len = 39;
        const real_t meanExcI_data[meanExcI_len] = { 8.3523e-05, 1.1718e-04, 6.4775e-05, 2.1155e-04, 2.6243e-04, 1.2896e-04, 1.8121e-04, 
                2.6380e-04, 4.1918e-04, 9.5147e-04, 0.0011, 2.6849e-04, 3.2329e-04, 3.8532e-04, 4.6027e-04, 5.5342e-04, 
                6.9002e-04, 9.2955e-04, 0.0014, 0.0028, 0.0029, 3.6888e-04, 4.2935e-04, 4.9667e-04, 5.7417e-04, 6.6360e-04, 
                7.7202e-04, 9.0685e-04, 0.0011, 0.0014, 0.0016, 0.0017, 0.0019, 0.0022, 0.0027, 0.0035, 0.0049, 0.0092, 0.0095};
        const len_t meanExcI_Zs[meanExcI_len] = { 2, 2, 3, 3, 3, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 10, 
                10, 10, 10, 10, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18};
        const len_t meanExcI_Z0s[meanExcI_len] = { 0, 1, 0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 
                4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};

        

        struct collqtyhand_settings *settings;

        gsl_integration_fixed_workspace **gsl_w = nullptr;

        virtual void InitializeGSLWorkspace();
    public:

        CollisionQuantityHandler(struct collqtyhand_settings *cq=nullptr);
        ~CollisionQuantityHandler();


        
        
        /**
         * The following three methods calculate and store all collision-frequency related quantities
         */
        virtual void CalculateCoulombLogarithms();            // lnL
        virtual void CalculateHiGiFuncs();                    // h_i and g_i
        virtual void CalculateCollisionFrequenciesFromHiGi(); // nu_s and nu_D
        
        virtual void CalculateIonisationRates();      // I, R and CE
        virtual void CalculateDerivedQuantities();    // Ec, Gamma_ava
        
        virtual real_t evaluatePsi0(len_t ir, real_t p);
        virtual real_t evaluatePsi1(len_t ir, real_t p);
        static real_t psi0Integrand(real_t s, void *params);
        static real_t psi1Integrand(real_t s, void *params);
        virtual real_t evaluateExp1OverThetaK(real_t Theta, real_t n);
        virtual real_t GetIonEffectiveSizeAj(len_t Z, len_t Z0);   // search atomic-data table for the Z, Z0 value. 
        virtual real_t GetMeanExcitationEnergy(len_t Z, len_t Z0); // search atomic-data table for the Z, Z0 value

        virtual void DeallocateLnLambdas();
        virtual void DeallocateHiGi();
        virtual void DeallocateCollisionFrequencies();
        virtual void DeallocateIonisationRates();
        virtual void DeallocateDerivedQuantities();
        virtual void DeallocateGSL();


        void SetEqSys(EquationSystem *es){
            this->eqSys = es;
        }
        virtual void Rebuild();

        void SetGrid(FVM::Grid *g, enum SimulationGenerator::momentumgrid_type mgtype){
            this->grid = g;
            this->gridtype = mgtype;
        }


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

        
         
        /**
         * NOTE: The below methods are not used in the standard DREAM workflow
         */

        
        //Calculate and stores nu_s and nu_D without storing h_i, g_i, lnL 
        virtual void CalculateCollisionFrequencies(); 
                

        
        
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


        




        virtual void SetIonSpecies(real_t **dens, len_t **Z, len_t **Z0, real_t *T);
        virtual void DeallocateIonSpecies();

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


        real_t *const* GetNuPar() const 
                { return this->collisionFrequencyNuPar; }
        const real_t  *GetNuPar(const len_t i) const 
                { return this->collisionFrequencyNuPar[i]; }
        real_t *const* GetNuPar_f1() const 
                { return this->collisionFrequencyNuPar_f1; }
        const real_t  *GetNuPar_f1(const len_t i) const 
                { return this->collisionFrequencyNuPar_f1[i]; }
        real_t *const* GetNuPar_f2() const 
                { return this->collisionFrequencyNuPar_f2; }
        const real_t  *GetNuPar_f2(const len_t i) const 
                { return this->collisionFrequencyNuPar_f2[i]; }


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
        
        
        real_t *const* GetHi(const len_t ir) const 
                { return this->HiFunc[ir]; }
        const real_t  GetHi(const len_t ir, const len_t i, const len_t j, const len_t iZ)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->HiFunc[ir][j*mg->GetNp1()+i][iZ]; }
        real_t *const* GetHi_f1(const len_t ir) const 
                { return this->HiFunc_f1[ir]; }
        const real_t  GetHi_f1(const len_t ir, const len_t i, const len_t j, const len_t iZ)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->HiFunc_f1[ir][j*mg->GetNp1()+i][iZ]; }
        real_t *const* GetHi_f2(const len_t ir) const 
                { return this->HiFunc_f2[ir]; }
        const real_t  GetHi_f2(const len_t ir, const len_t i, const len_t j, const len_t iZ)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->HiFunc_f2[ir][j*mg->GetNp1()+i][iZ]; }
        real_t *const* GetGi(const len_t ir) const 
                { return this->GiFunc[ir]; }
        const real_t  GetGi(const len_t ir, const len_t i, const len_t j, const len_t iZ)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->GiFunc[ir][j*mg->GetNp1()+i][iZ]; }
        real_t *const* GetGi_f1(const len_t ir) const 
                { return this->GiFunc_f1[ir]; }
        const real_t  GetGi_f1(const len_t ir, const len_t i, const len_t j, const len_t iZ)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->GiFunc_f1[ir][j*mg->GetNp1()+i][iZ]; }
        real_t *const* GetGi_f2(const len_t ir) const 
                { return this->GiFunc_f2[ir]; }
        const real_t  GetGi_f2(const len_t ir, const len_t i, const len_t j, const len_t iZ)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->GiFunc_f2[ir][j*mg->GetNp1()+i][iZ]; }
        const real_t *GetHCold(const len_t ir) const 
                { return this->HCold[ir]; }
        const real_t  GetHCold(const len_t ir, const len_t i, const len_t j)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->HCold[ir][j*mg->GetNp1()+i]; }
        const real_t *GetHCold_f1(const len_t ir) const 
                { return this->HCold_f1[ir]; }
        const real_t  GetHCold_f1(const len_t ir, const len_t i, const len_t j)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->HCold_f1[ir][j*mg->GetNp1()+i]; }
        const real_t *GetHCold_f2(const len_t ir) const 
                { return this->HCold_f2[ir]; }
        const real_t  GetHCold_f2(const len_t ir, const len_t i, const len_t j)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->HCold_f2[ir][j*mg->GetNp1()+i]; }
        const real_t *GetGCold(const len_t ir) const 
                { return this->GCold[ir]; }
        const real_t  GetGCold(const len_t ir, const len_t i, const len_t j)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->GCold[ir][j*mg->GetNp1()+i]; }
        const real_t *GetGCold_f1(const len_t ir) const 
                { return this->GCold_f1[ir]; }
        const real_t  GetGCold_f1(const len_t ir, const len_t i, const len_t j)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->GCold_f1[ir][j*mg->GetNp1()+i]; }
        const real_t *GetGCold_f2(const len_t ir) const 
                { return this->GCold_f2[ir]; }
        const real_t  GetGCold_f2(const len_t ir, const len_t i, const len_t j)  
                { FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir); return this->GCold_f2[ir][j*mg->GetNp1()+i]; }
        
    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP*/

    


