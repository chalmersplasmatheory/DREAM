#ifndef _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP
#define _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP

namespace DREAM { class EffectiveCriticalField; }



#include "DREAM/Equations/REPitchDistributionAveragedBACoeff.hpp"
#include "DREAM/Equations/AnalyticDistributionRE.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/IonHandler.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//#include "FVM/config.h"

namespace DREAM {
    class EffectiveCriticalField {
    public:
        

        struct ParametersForEceff {
            FVM::RadialGrid *rGrid; 
            SlowingDownFrequency *nuS; 
            PitchScatterFrequency *nuD; 
            FVM::fluxGridType fgType; 
            gsl_min_fminimizer *fmin; 
            CollisionQuantity::collqty_settings *collSettingsForEc; 
            OptionConstants::collqty_Eceff_mode Eceff_mode;
            IonHandler *ions;
            CoulombLogarithm *lnLambda;
            real_t thresholdToNeglectTrappedContribution;
        };


        /**
        * Parameter struct which is passed to all GSL functions involved in the Eceff calculations.
        */
        struct UContributionParams {
            SlowingDownFrequency *nuS;
            len_t ir;
            real_t Eterm; 
            gsl_min_fminimizer *fmin;
            real_t p_ex_lo;
            real_t p_ex_up; 
            real_t p_optimum;
            CollisionQuantity::collqty_settings *collSettingsForEc;
            AnalyticDistributionRE *analyticDist;
            REPitchDistributionAveragedBACoeff *EFieldTerm;
            REPitchDistributionAveragedBACoeff *SynchrotronTerm;
        };
        
    private:
        OptionConstants::collqty_Eceff_mode Eceff_mode;
        CollisionQuantity::collqty_settings *collSettingsForEc;

        FVM::RadialGrid *rGrid;
        len_t nr; 
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        IonHandler *ions;
        CoulombLogarithm *lnLambda;
        real_t thresholdToNeglectTrappedContribution;

        REPitchDistributionAveragedBACoeff *AveragedEFieldTerm;
        REPitchDistributionAveragedBACoeff *AveragedSynchrotronTerm;
        

        gsl_root_fdfsolver *fdfsolve;
        UContributionParams gsl_parameters;

        const real_t ECEFFOVERECTOT_INITGUESS = 1.0;
        const real_t POPTIMUM_INITGUESS = 10.0;
        
        real_t *ECRIT_ECEFFOVERECTOT_PREV = nullptr; // Eceff / Ectot in previous time step, used to accelerate Eceff algorithm
        real_t *ECRIT_POPTIMUM_PREV=nullptr;         // value of p which minimizes -U(p,Eceff)

    public:
        EffectiveCriticalField(ParametersForEceff*, AnalyticDistributionRE*);
        ~EffectiveCriticalField();

        bool GridRebuilt();
        void CalculateEffectiveCriticalField(const real_t *Ec_tot, const real_t *Ec_free, real_t *effectiveCriticalField);
        real_t CalculateEceffPPCFPaper(len_t ir);

        static real_t FindUExtremumAtE(real_t Eterm, void *par);
        static real_t FindUExtremumAtE_df(real_t Eterm, void *par);
        static void FindUExtremumAtE_fdf(real_t Eterm, void *par, real_t *f, real_t *df);
        static void FindPExInterval(
            real_t &p_ex_guess, real_t &p_ex_lower, real_t &p_ex_upper, 
            real_t &F_ex_guess, real_t &F_ex_lower, real_t &F_ex_upper,
            real_t p_upper_threshold, UContributionParams *params);
        static real_t UAtPFunc(real_t p, void *par); 
        void DeallocateQuantities();
    };
}

#endif/*_DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP*/