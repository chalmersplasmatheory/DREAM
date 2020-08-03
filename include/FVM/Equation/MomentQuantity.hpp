#ifndef _DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP
#define _DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


namespace DREAM::FVM {
    class MomentQuantity : public EquationTerm {
    protected:
        real_t *integrand;
        len_t nIntegrand = 0;

        FVM::Grid *fGrid;
        len_t momentId, fId;

        len_t nnz_per_row;
        void ResetIntegrand();
    public:
        MomentQuantity(Grid*, Grid*, len_t, len_t);
        virtual ~MomentQuantity();

        virtual len_t GetNumberOfNonZerosPerRow() const { return this->nnz_per_row; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const 
            { return GetNumberOfNonZerosPerRow(); }

        virtual bool GridRebuilt() override;

        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP*/
