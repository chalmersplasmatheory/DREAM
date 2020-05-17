#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP

namespace DREAM { class CollisionQuantityHandler; }

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"

#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/ParallelDiffusionFrequency.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"

#include <gsl/gsl_math.h>
#include "gsl/gsl_spline.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_interp2d.h>
#include <string>

namespace DREAM {
    class CollisionQuantityHandler{
    public:
        struct UExtremumParams {len_t ir; real_t Eterm; gsl_integration_workspace *gsl_w; CollisionQuantityHandler *collQtyHand;};
        void CalculatePStar();
        real_t evaluateBarNuSNuDAtP(len_t ir, real_t p){real_t p2=p*p; 
                return nuS->evaluateAtP(ir,p)*nuD->evaluateAtP(ir,p)*p2*p2*p2/(sqrt(1+p2)*(1+p2));}
        struct pStarFuncParams {real_t constTerm; len_t ir; CollisionQuantityHandler *collQtyHand;};
        static real_t pStarFunction(real_t, void *);
        void FindPInterval(len_t ir, real_t *p_lower, real_t *p_upper, pStarFuncParams pStar_params);
        void FindPStarRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func);

        void CalculateEffectiveCriticalField();
        real_t evaluateUAtP(len_t ir,real_t p, real_t Eterm,gsl_integration_workspace *gsl_ad_w);
        void FindECritInterval(len_t ir, real_t *E_lower, real_t *E_upper, UExtremumParams params);

        real_t evaluateComptonTotalCrossSectionAtP(real_t Eg, real_t pc);
        real_t evaluateComptonPhotonFluxSpectrum(real_t Eg);
        real_t evaluateBremsStoppingForceAtP(len_t ir,real_t p);
        

    private:
        const real_t constPreFactor = 4*M_PI
                                *Constants::r0*Constants::r0
                                *Constants::c;
        len_t nr;  // number of radial grid points 
        len_t nZ;  // number of atomic species
        len_t nzs; // number of ion species (including charge states)
        FVM::Grid *grid;
        //EquationSystem *eqSys = nullptr;
        FVM::UnknownQuantityHandler *unknowns = nullptr;
        IonHandler *ionHandler = nullptr;
        enum OptionConstants::momentumgrid_type gridtype;

        CoulombLogarithm *lnLambdaEE;
        CoulombLogarithm *lnLambdaEI;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        ParallelDiffusionFrequency *nuPar;


        // Ion densities on n x nZ
        real_t  *n_cold = nullptr;       // thermal free electron density in m^-3
        real_t  *T_cold = nullptr;                 // thermal free electron temperature in eV
        real_t *n_tot   = nullptr;
        real_t *E_term   = nullptr;
        const len_t  *ZAtomicNumber;          // atomic number (nuclear charge) of ion
        
        // Kinetic derived quantities on n 
        real_t *Ec_free=nullptr;        // Connor-Hastie field with only bound
        real_t *Ec_tot;                 // Connor-Hastie field with free+bound
        real_t *EDreic;                 // Dreicer field
        real_t *criticalREMomentum=nullptr; // Critical momentum for runaway p_star 
        real_t *avalancheRate;          // (dnRE/dt)_ava = nRE*Gamma_ava
        real_t *tritiumRate;            // (dnRE/dt)_Tritium
        real_t *comptonRate;            // (dnRE/dt)_Compton
        real_t *effectiveCriticalField; // Eceff: Gamma_ava(Eceff) = 0

        static const len_t  conductivityLenT;
        static const len_t  conductivityLenZ;
        static const real_t conductivityBraams[];
        static const real_t conductivityTmc2[];   // list of T/mc2 
        static const real_t conductivityX[];      // where X = 1/(1+Zeff) 
        

        struct CollisionQuantity::collqty_settings *settings;

        gsl_integration_fixed_workspace **gsl_w = nullptr;
        //gsl_interp_accel *gsl_acc  = gsl_interp_accel_alloc();

        const gsl_interp2d_type *gsl_T = gsl_interp2d_bilinear; 
        gsl_interp2d *gsl_cond = gsl_interp2d_alloc(gsl_T, 14,6);

        gsl_interp_accel *gsl_xacc = gsl_interp_accel_alloc();
        gsl_interp_accel *gsl_yacc = gsl_interp_accel_alloc();

        
        void CalculateDerivedQuantities();    // Ec, Gamma_ava

        void InitializeGSLWorkspace();
        real_t evaluateElectricalConductivity(len_t i);

        
        void CalculateGrowthRates();
        real_t evaluateTritiumRate(len_t ir);
        real_t evaluateComptonRate(len_t ir, gsl_integration_workspace *gsl_ad_w);
        void DeallocateDerivedQuantities();
    public:

        CollisionQuantityHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantity::collqty_settings *cqset);
        ~CollisionQuantityHandler();

        void gridRebuilt();

        void Rebuild();

        SlowingDownFrequency* GetNuS(){return nuS;}
        PitchScatterFrequency* GetNuD(){return nuD;}
        ParallelDiffusionFrequency* GetNuPar(){return nuPar;}
        CoulombLogarithm* GetLnLambdaEE(){return lnLambdaEE;}
        CoulombLogarithm* GetLnLambdaEI(){return lnLambdaEI;}

        const real_t GetLnLambdaC(len_t ir){return lnLambdaEE->GetLnLambdaC(ir);}
        const real_t *GetLnLambdaC(){return lnLambdaEE->GetLnLambdaC();}
        const real_t GetLnLambdaT(len_t ir){return lnLambdaEE->GetLnLambdaT(ir);}
        const real_t *GetLnLambdaT(){return lnLambdaEE->GetLnLambdaT();}
    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP*/

    


