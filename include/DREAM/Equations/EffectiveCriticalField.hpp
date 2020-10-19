#ifndef _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP
#define _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP

namespace DREAM { class EffectiveCriticalField; }

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/IonHandler.hpp"


//#include "FVM/config.h"

namespace DREAM {
    class EffectiveCriticalField {
    public:
        struct Param1 {FVM::RadialGrid *rGrid; SlowingDownFrequency *nuS; PitchScatterFrequency *nuD; FVM::fluxGridType fgType; 
                            gsl_integration_workspace *gsl_ad_w;
                            gsl_min_fminimizer *fmin; CollisionQuantity::collqty_settings *collSettingsForEc;};

        struct Param2 {CollisionQuantity::collqty_settings *collQtySettings;
                IonHandler *ions; gsl_root_fsolver *fsolve; OptionConstants::collqty_Eceff_mode Eceff_mode;};

        struct UContributionParams {FVM::RadialGrid *rGrid; SlowingDownFrequency *nuS; PitchScatterFrequency *nuD; len_t ir; real_t p; FVM::fluxGridType fgType; 
                            real_t Eterm; std::function<real_t(real_t,real_t,real_t)> Func; gsl_integration_workspace *gsl_ad_w;
                            gsl_min_fminimizer *fmin;real_t p_ex_lo; real_t p_ex_up; CollisionQuantity::collqty_settings *collSettingsForEc; int QAG_KEY;EffectiveCriticalField *ecEff;};
        
    private:
        // @@@Linnea temporarily put all parameters here! So we can have something that works
        len_t nr;
        CollisionQuantity::collqty_settings *collQtySettings;

        IonHandler *ions;

        gsl_root_fsolver *fsolve;
        OptionConstants::collqty_Eceff_mode Eceff_mode;

        FVM::RadialGrid *rGrid;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;

        FVM::fluxGridType fgType; // @@Linnea remove when working?
        gsl_integration_workspace *gsl_ad_w;
        
        gsl_min_fminimizer *fmin;
        
        CollisionQuantity::collqty_settings *collSettingsForEc;
        int QAG_KEY = GSL_INTEG_GAUSS31;

        
        //@@Linnea move def?

    public:
        EffectiveCriticalField(Param1*,Param2*);

        void CalculateEffectiveCriticalField(const real_t *Ec_tot, const real_t *Ec_free, real_t *effectiveCriticalField);

        static real_t FindUExtremumAtE(real_t Eterm, void *par);
        static void FindPExInterval(real_t *p_ex_guess, real_t *p_ex_lower, real_t *p_ex_upper, void *par, real_t p_upper_threshold);
        static real_t UAtPFunc(real_t p, void *par);

        //@@Linnea should be sent over to the runaway-fluid class for backwards compability.
        real_t testEvalU(len_t ir, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings);
        real_t evaluateAnalyticPitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm,CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w);
        real_t evaluateApproximatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm,CollisionQuantity::collqty_settings *inSettings);
        real_t evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w);

    };
}

#endif/*_DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP*/