#ifndef _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP
#define _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP

namespace DREAM { class EffectiveCriticalField; }


#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/AnalyticDistributionRE.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/IonHandler.hpp"


//#include "FVM/config.h"

namespace DREAM {
    class EffectiveCriticalField {
    public:
        struct ParametersForEceff {
            FVM::RadialGrid *rGrid; 
            SlowingDownFrequency *nuS; 
            PitchScatterFrequency *nuD; 
            FVM::fluxGridType fgType; 
            gsl_integration_workspace *gsl_ad_w;
            gsl_min_fminimizer *fmin; 
            CollisionQuantity::collqty_settings *collSettingsForEc, *collQtySettings; 
            gsl_root_fsolver *fsolve; 
            OptionConstants::collqty_Eceff_mode Eceff_mode;
            IonHandler *ions;
            CoulombLogarithm *lnLambda;
        };


        /**
        * Parameter struct which is passed to all GSL functions involved in the Eceff calculations.
        */
        struct UContributionParams {
            FVM::RadialGrid *rGrid;
            SlowingDownFrequency *nuS; PitchScatterFrequency *nuD;
            len_t ir;
            real_t p; 
            FVM::fluxGridType fgType;
            real_t Eterm;
            std::function<real_t(real_t,real_t,real_t)> Func; 
            gsl_integration_workspace *gsl_ad_w;
            gsl_min_fminimizer *fmin;
            real_t p_ex_lo;
            real_t p_ex_up; 
            real_t p_optimum;
            CollisionQuantity::collqty_settings *collSettingsForEc;
            int QAG_KEY;
            AnalyticDistributionRE *analyticDist;
        };
        
    private:
        OptionConstants::collqty_Eceff_mode Eceff_mode;
        CollisionQuantity::collqty_settings *collSettingsForEc;
        CollisionQuantity::collqty_settings *collQtySettings;

        FVM::RadialGrid *rGrid;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        IonHandler *ions;
        CoulombLogarithm *lnLambda;

        gsl_root_fsolver *fsolve;
        UContributionParams gsl_parameters;

        real_t *ECRIT_ECEFFOVERECTOT_PREV = nullptr; // Eceff / Ectot in previous time step, used to accelerate Eceff algorithm
        real_t *ECRIT_POPTIMUM_PREV=nullptr;         // value of p which minimizes -U(p,Eceff)

    public:
        EffectiveCriticalField(ParametersForEceff*, AnalyticDistributionRE*);
        ~EffectiveCriticalField();

        bool GridRebuilt();
        void CalculateEffectiveCriticalField(const real_t *Ec_tot, const real_t *Ec_free, real_t *effectiveCriticalField);
        real_t CalculateEceffPPCFPaper(len_t ir);

        static real_t FindUExtremumAtE(real_t Eterm, void *par);
        static void FindPExInterval(real_t *p_ex_guess, real_t *p_ex_lower, real_t *p_ex_upper, real_t p_upper_threshold, 
        UContributionParams *params);
        static real_t UAtPFunc(real_t p, void *par);
    };
}

#endif/*_DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP*/