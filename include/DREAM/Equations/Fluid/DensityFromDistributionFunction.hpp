#ifndef _DREAM_EQUATION_FLUID_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/MomentQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class DensityFromDistributionFunction : public FVM::MomentQuantity {
    public:
        // Specifies unit in which pThreshold is given to constructor
        enum pThresholdMode {
            P_THRESHOLD_MODE_MC,
            P_THRESHOLD_MODE_THERMAL
        };
    
    private:
        real_t pThreshold;
        pThresholdMode pMode;

        real_t p0;

        real_t GetThreshold(len_t ir, FVM::UnknownQuantityHandler*);
    public:
        DensityFromDistributionFunction(
            FVM::Grid*, FVM::Grid*, len_t, len_t,
            real_t pThreshold = 0, pThresholdMode pMode = P_THRESHOLD_MODE_MC
        );
        virtual ~DensityFromDistributionFunction();

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP*/
