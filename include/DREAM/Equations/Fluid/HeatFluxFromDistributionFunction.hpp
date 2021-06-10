#ifndef _DREAM_EQUATION_FLUID_HEAT_FLUX_FROM_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_HEAT_FLUX_FROM_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/MomentQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM {
    class HeatFluxFromDistributionFunction : public FVM::MomentQuantity {
    private:
        real_t scaleFactor;
    public:
        HeatFluxFromDistributionFunction(
            FVM::Grid*, FVM::Grid*, len_t, len_t, FVM::UnknownQuantityHandler*,
            real_t pThreshold = 0, pThresholdMode pMode = FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_MC,
            real_t scaleFactor = 1.0);
        virtual ~HeatFluxFromDistributionFunction(){}

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
    };
}

#endif/*_DREAM_EQUATION_FLUID_HEAT_FLUX_FROM_DISTRIBUTION_FUNCTION_HPP*/
