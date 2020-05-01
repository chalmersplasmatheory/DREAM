#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP

namespace DREAM { class CollisionQuantityHandler; }

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"
#include <gsl/gsl_math.h>
#include "gsl/gsl_spline.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_interp2d.h>
#include <string>

namespace DREAM {
    class CollisionQuantityHandler{

    public:
        struct collqtyhand_settings {
            enum OptionConstants::collqty_collfreq_type 
                        collfreq_type   = OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED;
            enum OptionConstants::collqty_collfreq_mode 
                        collfreq_mode   = OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
            enum OptionConstants::collqty_lnLambda_type 
                        lnL_type        = OptionConstants::COLLQTY_LNLAMBDA_CONSTANT;
            enum OptionConstants::uqty_n_cold_eqn       
                        ncold_type      = OptionConstants::UQTY_N_COLD_EQN_PRESCRIBED;
            enum OptionConstants::eqterm_nonlinear_mode
                        nonlinear_mode  = OptionConstants::EQTERM_NONLINEAR_MODE_NEGLECT;
        };

    private:
        const real_t constPreFactor = 4*M_PI
                                *Constants::r0*Constants::r0
                                *Constants::c;
        len_t n;   // number of "radial grid points" (or sets of ion species) 
        len_t nZ;  // number of atomic species
        len_t nzs; // number of ion species (including charge states)
        FVM::Grid *grid;
        //EquationSystem *eqSys = nullptr;
        FVM::UnknownQuantityHandler *unknowns = nullptr;
        IonHandler *ionHandler = nullptr;
        enum OptionConstants::momentumgrid_type gridtype;

        // Ion densities on n x nZ
        real_t  *n_cold = nullptr;       // thermal free electron density in m^-3
        real_t  *T_cold;                 // thermal free electron temperature in eV
        //real_t **ionDensity=nullptr;     // ion densities in m^-3
        const len_t  *ZAtomicNumber;          // atomic number (nuclear charge) of ion
        //const len_t  *Z0ChargeNumber;         // charge number (net charge) of ion
        
        real_t *n_tot = nullptr;

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

/*
        // Ionisation and recombination on n x nZ
        real_t **ionisationRateCold=nullptr;  // ionisation rate by thermal cold electrons
        real_t **ionisationRateHot;       // ionisation rate by kinetic electrons
        //         ^-- perhaps this should be a function of p? (or even f?)
        //             unreasonable to store on n x nZ x (np1 x np2)..?
        real_t **ionisationRateREFluid;            // ionisation rate by RE fluid
        real_t **recombinationRateRadiative;  // radiative recombination rate
        real_t **chargeExchangeZP;            // impurity-proton charge exchange rates
        real_t  *chargeExchangeHP;            // hydrogen-proton charge exchange rate
*/        
        
        // Atomic parameters on nZ
        real_t *ionisationPotential=nullptr;  // Ionisation energies loaded from file
        real_t *meanExcitationEnergy;         // For nu_s. Loaded from file or calculated.
        real_t *ionEffectiveSizeAj;           // For nu_D. Loaded from file or calculated.

        // Kinetic derived quantities on n 
        real_t *Ec_free=nullptr;        // Connor-Hastie field with only bound
        real_t *Ec_tot;                 // Connor-Hastie field with free+bound
        real_t *EDreic;                 // Dreicer field
        real_t *criticalREMomentum=nullptr; // Critical momentum for runaway p_star 
        real_t *avalancheRate;          // (dnRE/dt)_ava = nRE*Gamma_ava
        real_t *tritiumRate;            // (dnRE/dt)_Tritium
        real_t *comptonRate;            // (dnRE/dt)_Compton
        real_t *effectiveCriticalField; // Eceff: Gamma_ava(Eceff) = 0

        
        // atomic data in no particular order, but ...data[i] corresponds to the value for charge Z = ...Zs[i] and Z0 = ...Z0s[i]
        static const len_t  ionSizeAj_len; 
        static const real_t ionSizeAj_data[];       
        static const real_t ionSizeAj_Zs[];
        static const real_t ionSizeAj_Z0s[];

        static const len_t  meanExcI_len;
        static const real_t meanExcI_data[];
        static const real_t meanExcI_Zs[];
        static const real_t meanExcI_Z0s[];
        
        static const len_t  conductivityLenT;
        static const len_t  conductivityLenZ;
        static const real_t conductivityBraams[];
        static const real_t conductivityTmc2[];   // list of T/mc2 
        static const real_t conductivityX[];      // where X = 1/(1+Zeff) 
        

        real_t **nonlinearApMat = nullptr;
        real_t **nonlinearDppMat;
        real_t **nonlinearNuDMat;


        struct collqtyhand_settings *settings;

        gsl_integration_fixed_workspace **gsl_w = nullptr;
        //gsl_interp_accel *gsl_acc  = gsl_interp_accel_alloc();

        const gsl_interp2d_type *gsl_T = gsl_interp2d_bilinear; 
        gsl_interp2d *gsl_cond = gsl_interp2d_alloc(gsl_T, 14,6);

        gsl_interp_accel *gsl_xacc = gsl_interp_accel_alloc();
        gsl_interp_accel *gsl_yacc = gsl_interp_accel_alloc();

        virtual void InitializeGSLWorkspace();

                
        static const len_t numSupportedSpecies;
        static const len_t Zdata[];
        static const std::string stringsdata[];
        //virtual void ReadADASDataFromFile(std::string coefficientType, len_t Z, real_t *&log10n, real_t *&log10T, real_t **&Coefficients);
        //virtual std::string GetADASPath(std::string coefficientType, len_t Z);
    public:

        CollisionQuantityHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct collqtyhand_settings *cqset);
        virtual ~CollisionQuantityHandler();


        
        
        /**
         * The following three methods calculate and store all collision-frequency related quantities
         */
        virtual void CalculateCoulombLogarithms();            // lnL
        virtual void CalculateHiGiFuncs();                    // h_i and g_i
        virtual void CalculateCollisionFrequenciesFromHiGi(); // nu_s and nu_D
        
        //virtual void CalculateIonisationRates();      // I, R and CE
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
        //virtual void DeallocateIonisationRates();
        virtual void DeallocateDerivedQuantities();
        virtual void DeallocateGSL();
        virtual void DeallocateNonlinearMatrices();

        virtual void Rebuild();
/*
        void SetUnknowns(FVM::UnknownQuantityHandler *u){
            this->unknowns = u;
            // this->nZ = number of ion species stored
        }
        

        void SetGrid(FVM::Grid *g, enum OptionConstants::momentumgrid_type mgtype){
            this->grid = g;
            this->gridtype = mgtype;
            this->n = g->GetNr();
        }
*/

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
        
        virtual real_t evaluateLnLambdaEEAtP(len_t i,real_t p);
        virtual real_t evaluateLnLambdaEIAtP(len_t i,real_t p);
        virtual real_t evaluateLnLambdaC(len_t i);

        virtual real_t evaluateElectricalConductivity(len_t i);

        virtual void CalculatePStar();
        virtual real_t evaluateBarNuSNuDAtP(len_t ir, real_t p)
                {return evaluateNuSAtP(ir,p)*evaluateNuDAtP(ir,p) *p*p*p*p*p*p/(sqrt(1+p*p)*(1+p*p));}
        struct pStarFuncParams {real_t constTerm; len_t ir; CollisionQuantityHandler *collQtyHand;};
        static real_t pStarFunction(real_t, void *);
        virtual void FindPInterval(len_t ir, real_t *p_lower, real_t *p_upper, pStarFuncParams pStar_params);
        virtual void FindPStarRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func);

        virtual void CalculateEffectiveCriticalField();
        virtual real_t evaluateUAtP(len_t ir,real_t p, real_t Eterm,gsl_integration_workspace *gsl_ad_w);
        struct UExtremumParams {len_t ir; real_t Eterm; gsl_integration_workspace *gsl_w; CollisionQuantityHandler *collQtyHand;};
        virtual void FindECritInterval(len_t ir, real_t *E_lower, real_t *E_upper, UExtremumParams params);
        
        virtual void CalculateGrowthRates();
        virtual real_t evaluateTritiumRate(len_t ir);
        virtual real_t evaluateComptonRate(len_t ir, gsl_integration_workspace *gsl_ad_w);
        virtual real_t evaluateComptonTotalCrossSectionAtP(real_t Eg, real_t pc);
        virtual real_t evaluateComptonPhotonFluxSpectrum(real_t Eg);


        virtual void calculateIsotropicNonlinearOperatorMatrices();
        virtual void addNonlinearContribNuS(real_t**&);
        virtual void addNonlinearContribNuPar(real_t**&);
        virtual void addNonlinearContribNuD(real_t**&);
        
        /**
         * NOTE: The below methods are not used in the standard DREAM workflow
         */

        
        //Calculate and stores nu_s and nu_D without storing h_i, g_i, lnL 
        virtual void CalculateCollisionFrequencies(); 
                



        




        //virtual void SetIonSpecies(real_t **dens, len_t *Z, len_t *Z0, real_t *T);
        //virtual void DeallocateIonSpecies();

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


        
        //virtual void LoadAtomicData();

        void SetAtomicParameters(real_t *I, real_t *Imean){
            this->meanExcitationEnergy = Imean;
            this->ionisationPotential = I;
        }

/*
        void SetIonisationRates(real_t **Icold, real_t **Ihot, real_t **IRE,
                                    real_t **RR, real_t **CEZP, real_t *CEHP){
            DeallocateIonisationRates();
            this->ionisationRateCold = Icold;
            this->ionisationRateHot = Ihot;
            this->ionisationRateREFluid = IRE;
            this->recombinationRateRadiative = RR;
            this->chargeExchangeZP = CEZP;
            this->chargeExchangeHP = CEHP;       
        }
*/

        void SetDerivedQuantities(real_t *Ec, real_t *Ectot, real_t *ED, 
                                real_t *Gamma,  real_t *Eceff, real_t *pStar){ 
            this->Ec_free = Ec;
            this->Ec_tot  = Ectot;
            this->EDreic  = ED;
            this->criticalREMomentum  = pStar;
            this->avalancheRate = Gamma;
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

    


