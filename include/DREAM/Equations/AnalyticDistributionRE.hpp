#ifndef _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP
#define _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP

namespace DREAM { class AnalyticDistributionRE; }

#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/IonHandler.hpp"


//#include "FVM/config.h"

namespace DREAM {
    class AnalyticDistributionRE {
        private:
            FVM::RadialGrid *rGrid;
            PitchScatterFrequency *nuD; 
            OptionConstants::collqty_Eceff_mode Eceff_mode;
        public:
            AnalyticDistributionRE(FVM::RadialGrid *rGrid, PitchScatterFrequency *nuD, OptionConstants::collqty_Eceff_mode Eceff_mode);
            real_t evaluateAnalyticPitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings, 
            gsl_integration_workspace *gsl_ad_w);
            real_t evaluateApproximatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm,CollisionQuantity::collqty_settings *inSettings);
            real_t evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings, 
            gsl_integration_workspace *gsl_ad_w);
    };
}

#endif/*_DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_RE_HPP*/