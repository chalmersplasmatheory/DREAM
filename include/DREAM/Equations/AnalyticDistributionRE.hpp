#ifndef _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP
#define _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP

namespace DREAM { class AnalyticDistributionRE; }

#include "DREAM/Equations/AnalyticDistribution.hpp"
#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

namespace DREAM {
    class AnalyticDistributionRE : public AnalyticDistribution {
        public:
            enum dist_mode {
                RE_PITCH_DIST_SIMPLE = 1,
                RE_PITCH_DIST_FULL = 2
            };
        private:
            PitchScatterFrequency *nuD; 
            CollisionQuantity::collqty_settings *collSettings;
            dist_mode mode;
            real_t thresholdToNeglectTrappedContribution;

            len_t id_Eterm;
            len_t id_nre;
            gsl_integration_workspace *gsl_ad_w;

            gsl_spline **xi0OverXiSpline = nullptr;
            gsl_interp_accel **xiSplineAcc = nullptr;
            static const len_t N_SPLINE = 10; 
            real_t *integralOverFullPassing = nullptr;
            
            gsl_interp_accel **REDistNormFactor_Accel = nullptr;
            gsl_spline **REDistNormFactor_Spline = nullptr;

            void Deallocate();
            void constructXiSpline();
            void constructVpSplines();
        public:
            AnalyticDistributionRE(
                FVM::RadialGrid*, FVM::UnknownQuantityHandler*, PitchScatterFrequency*, 
                CollisionQuantity::collqty_settings*, dist_mode, 
                real_t thresholdToNeglectTrappedContribution
            );
            virtual ~AnalyticDistributionRE();

            virtual bool GridRebuilt() override;

            real_t GetAatP(len_t ir, real_t p, CollisionQuantity::collqty_settings *set=nullptr, real_t *Eterm=nullptr);

            virtual real_t evaluateEnergyDistribution(len_t ir, real_t p, real_t *dfdp=nullptr, real_t *dfdr=nullptr) override;
            virtual real_t evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) override;

            real_t evaluateAnalyticPitchDistributionFromA(
                len_t ir, real_t xi0, real_t A
            );

            real_t EvaluateVpREAtA(len_t ir, real_t A);

            real_t evaluateApproximatePitchDistributionFromA(len_t ir, real_t xi0, real_t A);
            real_t evaluatePitchDistributionFromA(len_t ir, real_t xi0, real_t A);

            FVM::RadialGrid *GetRadialGrid() {return this->rGrid;}

    };
}

#endif/*_DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP*/
