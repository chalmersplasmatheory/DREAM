#ifndef _DREAM_EQUATIONS_FLUID_SOURCE_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_SOURCE_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
namespace DREAM {
    class FluidSourceTerm
        : public FVM::EquationTerm {

    private:
        real_t *sourceVec = nullptr;
    protected:
        FVM::UnknownQuantityHandler *unknowns;
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) = 0;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) = 0;
        virtual void NormalizeSourceToConstant(const real_t c, real_t *normFactors = nullptr);
    public:
        FluidSourceTerm(FVM::Grid*, FVM::UnknownQuantityHandler*);
        ~FluidSourceTerm();

        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x) override;

        virtual len_t GetNumberOfNonZerosPerRow() const override 
            { return 1; }

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual bool GridRebuilt() override;
    };
}

#endif/*_DREAM_EQUATIONS_FLUID_SOURCE_TERM_HPP*/


