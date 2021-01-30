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
            CollisionQuantity::collqty_settings *collSettingsForEc; 
            gsl_root_fdfsolver *fdfsolve; 
            OptionConstants::collqty_Eceff_mode Eceff_mode;
            IonHandler *ions;
            CoulombLogarithm *lnLambda;
            real_t thresholdToNeglectTrappedContribution;
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
            real_t A;
            real_t(*BAFunc)(real_t,real_t,real_t,real_t,void*);
            void *BAFunc_par;
            const int_t *BAList;
            std::function<real_t(real_t)> preFactorFunc; 
            gsl_integration_workspace *gsl_ad_w;
            gsl_integration_workspace *gsl_ad_w2;
            gsl_min_fminimizer *fmin;
            real_t p_ex_lo;
            real_t p_ex_up; 
            real_t p_optimum;
            CollisionQuantity::collqty_settings *collSettingsForEc;
            AnalyticDistributionRE *analyticDist;
            real_t CONST_E;
            real_t CONST_EFact;
            real_t CONST_Synch;
            gsl_spline **EContribSpline; 
            gsl_spline **SynchContribSpline;
            gsl_interp_accel **EContribAcc;
            gsl_interp_accel **SynchContribAcc;
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

        gsl_root_fdfsolver *fdfsolve;
        UContributionParams gsl_parameters;

        static const len_t N_A_VALUES = 15; 
        real_t X_vec[N_A_VALUES];
        real_t **EOverUnityContrib=nullptr;
        real_t **SynchOverUnityContrib=nullptr;
        real_t *ECRIT_ECEFFOVERECTOT_PREV = nullptr; // Eceff / Ectot in previous time step, used to accelerate Eceff algorithm
        real_t *ECRIT_POPTIMUM_PREV=nullptr;         // value of p which minimizes -U(p,Eceff)

        /** 
         * The splines are stored and evaluated in the
         * variable X = A/(1+A) which maps the interval
         * [0,inf] in A to [0,1] in X. These functions
         * convert between A and X.
         */
        static real_t GetAFromX(real_t X)
            { return sqrt(X)/(1.0-sqrt(X)); }
        static real_t GetXFromA(real_t A)
            { real_t x = A/(1+A); return x*x; }

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
        void CreateLookUpTableForUIntegrals(UContributionParams *par, real_t &EContrib, real_t &SynchContrib);
        void DeallocateQuantities();
    };
}

#endif/*_DREAM_EQUATIONS_EFFECTIVE_CRITICAL_FIELD_HPP*/