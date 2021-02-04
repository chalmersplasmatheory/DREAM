#ifndef _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_HOTTAIL_HPP
#define _DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_HOTTAIL_HPP

namespace DREAM { class AnalyticDistributionHottail; }

#include "DREAM/Equations/AnalyticDistribution.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
namespace DREAM {
    class AnalyticDistributionHottail : public AnalyticDistribution {
        private:
            OptionConstants::uqty_f_hot_dist_mode type;
            real_t *n0, *T0;
            len_t id_tau;

            real_t 
                *preFactor=nullptr,
                *betaTh = nullptr;

        public:
            AnalyticDistributionHottail(FVM::RadialGrid*, FVM::UnknownQuantityHandler*, real_t*, real_t*, OptionConstants::uqty_f_hot_dist_mode);
            virtual ~AnalyticDistributionHottail();
            virtual bool GridRebuilt() override;

            const real_t GetInitialThermalMomentum(const len_t ir) const {return betaTh[ir];}
            const real_t* GetInitialThermalMomentum() const {return betaTh;}

            virtual real_t evaluateEnergyDistribution(len_t ir, real_t p, real_t *dfdp=nullptr, real_t *dfdr=nullptr) override;            
            real_t evaluateEnergyDistributionFromTau(len_t ir, real_t p, real_t tau, real_t *dfdp=nullptr, real_t *dfdr=nullptr, real_t *dFdpOverF=nullptr, real_t *dFdTau=nullptr);

            // isotropic pitch distribution
            virtual real_t evaluatePitchDistribution(len_t /*ir*/, real_t /*xi0*/, real_t /*p*/, real_t *dfdxi0=nullptr, real_t *dfdp=nullptr, real_t *dfdr=nullptr) override {
                if(dfdxi0!=nullptr)
                    *dfdxi0 = 0;
                if(dfdp!=nullptr)
                    *dfdp = 0;
                if(dfdr!=nullptr)
                    *dfdr = 0;

                return 1.0/(4*M_PI);
            }
    };
}

#endif/*_DREAM_EQUATIONS_ANALYTIC_DISTRIBUTION_HOTTAIL_HPP*/
