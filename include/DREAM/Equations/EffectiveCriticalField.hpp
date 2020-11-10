#ifndef _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP
#define _DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP

namespace DREAM { class EffectiveCriticalField; }


#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/AnalyticDistributionRE.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/IonHandler.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h> // remove when cleaned up!
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
            real_t A;// remove p and Eterm later?? @@Linnea
            std::function<real_t(real_t,real_t,real_t)> Func; 
            gsl_integration_workspace *gsl_ad_w;
            gsl_min_fminimizer *fmin;
            real_t p_ex_lo;
            real_t p_ex_up; 
            CollisionQuantity::collqty_settings *collSettingsForEc;
            int QAG_KEY;
            AnalyticDistributionRE *analyticDist;
            real_t *A_vec;
            len_t N_A_VALUES;
            real_t **EContribIntegral;
            real_t **SynchContribIntegral;
            real_t **UnityContribIntegral;
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

        //real_t *pc_cylApprox;
        //real_t *Eceff_cylApprox;

        gsl_root_fsolver *fsolve;
        UContributionParams gsl_parameters;

        static const len_t N_A_VALUES = 50;
        real_t A_vec[N_A_VALUES];
        real_t **EContribIntegral=nullptr;
        real_t **SynchContribIntegral=nullptr;
        real_t **UnityContribIntegral=nullptr;

    public:
        EffectiveCriticalField(ParametersForEceff*, AnalyticDistributionRE*);
        ~EffectiveCriticalField();

        void CalculateEffectiveCriticalField(const real_t *Ec_tot, const real_t *Ec_free, real_t *effectiveCriticalField);
        real_t CalculateEceffPPCFPaper(len_t ir);

        static real_t FindUExtremumAtE(real_t Eterm, void *par);
        static void FindPExInterval(real_t *p_ex_guess, real_t *p_ex_lower, real_t *p_ex_upper, real_t p_upper_threshold, 
        UContributionParams *params);
        static real_t UAtPFuncNoSpline(real_t p, void *par);
        static real_t UAtPFunc(real_t p, void *par); // later: maybe it doesn't need to be static; different input too?
        void CreateLookUpTableForUIntegrals(UContributionParams *par, real_t *EContrib, real_t *UnityContrib, real_t *SynchContrib);
        static real_t InterpolateVector(real_t *x, real_t *y, real_t xNew, len_t Nx);
        // new stuff!
        static real_t GetU2atPandE(const gsl_vector *v, void *par);
    };
}

#endif/*_DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP*/