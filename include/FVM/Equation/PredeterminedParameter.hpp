#ifndef _DREAM_FVM_EQUATION_PREDETERMINED_PARAMETER_HPP
#define _DREAM_FVM_EQUATION_PREDETERMINED_PARAMETER_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class PredeterminedParameter : public EquationTerm {
    protected:
        real_t *currentData = nullptr;

    public:
        PredeterminedParameter(Grid *g);
        virtual ~PredeterminedParameter();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual const real_t *GetData() { return this->currentData; }
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_PREDETERMINED_PARAMETER_HPP*/
