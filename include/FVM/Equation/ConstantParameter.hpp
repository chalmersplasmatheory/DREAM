#ifndef _DREAM_FVM_EQUATION_NULL_PARAMETER_HPP
#define _DREAM_FVM_EQUATION_NULL_PARAMETER_HPP

#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class ConstantParameter : public PredeterminedParameter {
    public:
        ConstantParameter(Grid *g, const real_t v);
        virtual ~ConstantParameter();

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override {}
    };
}

#endif/*_DREAM_FVM_EQUATION_NULL_PARAMETER_HPP*/
