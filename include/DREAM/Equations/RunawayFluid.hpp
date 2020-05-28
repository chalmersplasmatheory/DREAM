#ifndef _DREAM_EQUATIONS_RUNAWAY_FLUID_HPP
#define _DREAM_EQUATIONS_RUNAWAY_FLUID_HPP

namespace DREAM { class RunawayFluid; }
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"

#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

namespace DREAM {
    class RunawayFluid{
    private:
        const real_t constPreFactor = 4*M_PI
                                *Constants::r0*Constants::r0
                                *Constants::c;

        FVM::RadialGrid *rGrid;
        FVM::UnknownQuantityHandler *unknowns;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        CoulombLogarithm *lnLambdaEE;
        len_t nr;

        len_t id_ncold;
        len_t id_ntot;
        len_t id_Tcold;
        len_t id_Eterm;

        real_t *ncold;
        real_t *ntot;
        real_t *Tcold;
        real_t *Eterm;

        real_t *Ec_free=nullptr;                // Connor-Hastie field with only bound
        real_t *Ec_tot=nullptr;                 // Connor-Hastie field with free+bound
        real_t *EDreic=nullptr;                 // Dreicer field
        real_t *criticalREMomentum=nullptr;     // Critical momentum for runaway p_star 
        real_t *pc_COMPLETESCREENING = nullptr;
        real_t *pc_NOSCREENING = nullptr;
        real_t *avalancheGrowthRate=nullptr;    // (dnRE/dt)_ava = nRE*Gamma_ava
        real_t *tritiumRate=nullptr;            // (dnRE/dt)_Tritium
        real_t *comptonRate=nullptr;            // (dnRE/dt)_Compton
        real_t *effectiveCriticalField=nullptr; // Eceff: Gamma_ava(Eceff) = 0

        //real_t *exactPitchDistNormalization;
        //real_t *approximatePitchDistNormalization;

       


        bool gridRebuilt;
        void AllocateQuantities();
        void DeallocateQuantities();
        //void AllocateGSLWorkspaces();
        //void FreeGSLWorkspaces();

        void CalculateDerivedQuantities();
        void CalculateEffectiveCriticalField(bool useApproximateMethod);
        void CalculateCriticalMomentum();

        static void FindECritInterval(len_t ir, real_t *E_lower, real_t *E_upper, void *par);
        static void FindPExInterval(real_t *p_ex_guess, real_t *p_ex_lower, real_t *p_ex_upper, void *par, real_t p_upper_threshold);
        static void FindRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func, gsl_root_fsolver *s);
//        void CalculateDistributionNormalizationFactors();

        real_t BounceAverageFunc(len_t ir, std::function<real_t(real_t,real_t)> Func);

        static real_t FindUExtremumAtE(real_t Eterm, void *par);
        static real_t evaluateNegUAtP(real_t p, void *par);
        static real_t evaluateApproximateUAtP(real_t p, void *par);
        static real_t UAtPFunc(real_t p, void *par);
//        real_t evaluateDistNormalization(len_t ir, real_t p, real_t Eterm);
//        real_t evaluateApproximateDistNormalization(len_t ir, real_t p, real_t Eterm);
    protected:
    public:
        RunawayFluid(FVM::Grid *g, FVM::UnknownQuantityHandler *u, SlowingDownFrequency *nuS, 
        PitchScatterFrequency *nuD, CoulombLogarithm *lnLambdaEE);
        ~RunawayFluid();

        real_t testEvalU(len_t ir, real_t p, real_t Eterm, bool useApproximateMethod);

        real_t evaluateAnalyticPitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, gsl_integration_workspace *gsl_ad_w);
        real_t evaluateApproximatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm);

        void Rebuild(bool useApproximateMethod);
        void GridRebuilt(){gridRebuilt = true;}
        const real_t GetEffectiveCriticalField(len_t ir) const
            {return effectiveCriticalField[ir];}
        const real_t* GetEffectiveCriticalField() const
            {return effectiveCriticalField;}
        
        const real_t GetDreicerElectricField(len_t ir) const
            {return EDreic[ir];}
        const real_t* GetDreicerElectricField() const
            {return EDreic;}
        
        const real_t GetConnorHastieField_COMPLETESCREENING(len_t ir) const
            {return Ec_free[ir];}
        const real_t* GetConnorHastieField_COMPLETESCREENING() const
            {return Ec_free;}
        
        const real_t GetConnorHastieField_NOSCREENING(len_t ir) const
            {return Ec_tot[ir];}
        const real_t* GetConnorHastieField_NOSCREENING() const
            {return Ec_tot;}
        
        

        const real_t GetAvalancheGrowthRate(len_t ir) const
            {return avalancheGrowthRate[ir];}
        const real_t* GetAvalancheGrowthRate() const
            {return avalancheGrowthRate;}
        
        const real_t GetTritiumRunawayRate(len_t ir) const
            {return tritiumRate[ir];}
        const real_t* GetTritiumRunawayRate() const
            {return tritiumRate;}
        
        const real_t GetComptonRunawayRate(len_t ir) const
            {return comptonRate[ir];}
        const real_t* GetComptonRunawayRate() const
            {return comptonRate;}

        const real_t GetEffectiveCriticalRunawayMomentum(len_t ir) const
            {return criticalREMomentum[ir];}
        const real_t* GetEffectiveCriticalRunawayMomentum() const
            {return criticalREMomentum;}
        

    };

}


#endif/*_DREAM_EQUATIONS_RUNAWAY_FLUID_HPP*/

    


