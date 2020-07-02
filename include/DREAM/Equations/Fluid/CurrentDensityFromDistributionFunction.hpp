#ifndef _DREAM_EQUATION_FLUID_CURRENT_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_CURRENT_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/MomentQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM {
    class CurrentDensityFromDistributionFunction : public FVM::MomentQuantity {
    public:
        CurrentDensityFromDistributionFunction(FVM::Grid*, FVM::Grid*, len_t, len_t);
        virtual ~CurrentDensityFromDistributionFunction();

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
    };
}

#endif/*_DREAM_EQUATION_FLUID_CURRENT_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP*/
