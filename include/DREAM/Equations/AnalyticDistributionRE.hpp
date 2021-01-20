#ifndef _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP
#define _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP

namespace DREAM { class AnalyticDistributionRE; }

#include "DREAM/Equations/AnalyticDistribution.hpp"
#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include <gsl/gsl_integration.h>

namespace DREAM {
    class AnalyticDistributionRE : public AnalyticDistribution {
        public:
            enum dist_mode {
                RE_PITCH_DIST_SIMPLE = 1,
                RE_PITCH_DIST_FULL = 2
            };
        private:
            FVM::RadialGrid *rGrid;
            PitchScatterFrequency *nuD; 
            CollisionQuantity::collqty_settings *collSettings;
            dist_mode mode;
            real_t thresholdToNeglectTrappedContribution;

            gsl_integration_workspace *gsl_ad_w;
            len_t id_Eterm;
        public:
            AnalyticDistributionRE(
                FVM::RadialGrid*, FVM::UnknownQuantityHandler*, PitchScatterFrequency*, 
                CollisionQuantity::collqty_settings*, dist_mode, 
                real_t thresholdToNeglectTrappedContribution
            );
            virtual ~AnalyticDistributionRE();

            virtual real_t evaluateFullDistribution(len_t ir, real_t xi0, real_t p, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) override;
            virtual real_t evaluateEnergyDistribution(len_t ir, real_t p, real_t *dfdp=nullptr, real_t *dfdr=nullptr) override;
            virtual real_t evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) override;

            real_t evaluateAnalyticPitchDistributionFromA(
                len_t ir, real_t xi0, real_t A
            );
            real_t evaluateApproximatePitchDistributionFromA(len_t ir, real_t xi0, real_t A);
            real_t evaluatePitchDistributionFromA(len_t ir, real_t xi0, real_t A);
    };
}

#endif/*_DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP*/
