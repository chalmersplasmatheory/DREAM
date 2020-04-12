#ifndef _DREAM_EQUATION_FLUID_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/MomentQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class DensityFromDistributionFunction : public FVM::MomentQuantity {
    public:
        DensityFromDistributionFunction(FVM::Grid*, FVM::Grid*, len_t, len_t);
        virtual ~DensityFromDistributionFunction();

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
    };
}

#endif/*_DREAM_EQUATION_FLUID_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP*/
