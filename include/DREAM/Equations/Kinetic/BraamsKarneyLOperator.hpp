#ifndef _DREAM_EQUATIONS_BRAAMSKARNEY_L_OPERATOR_HPP
#define _DREAM_EQUATIONS_BRAAMSKARNEY_L_OPERATOR_HPP

#include "FVM/Equation/EquationTerm.hpp"

namespace DREAM {

    class BraamsKarneyLOperator : public FVM::EquationTerm {
    private:
        real_t a;

        template<typename F1>
        void SetElements(F1 X);

    public:
        BraamsKarneyLOperator(FVM::Grid *grid, real_t a);

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
        virtual len_t GetNumberOfNonZerosPerRow() const override;

        bool SetJacobianBlock(const len_t unknId, const len_t derivId, FVM::Matrix *jac, const real_t *) override;
        void SetMatrixElements(FVM::Matrix *mat, real_t*) override;
        void SetVectorElements(real_t *vec, const real_t *f) override;
    };

}

#endif/*_DREAM_EQUATIONS_BRAAMSKARNEY_RHS_HPP*/
