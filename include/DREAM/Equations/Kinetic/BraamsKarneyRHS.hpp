#ifndef _DREAM_EQUATIONS_BRAAMSKARNEY_RHS_HPP
#define _DREAM_EQUATIONS_BRAAMSKARNEY_RHS_HPP

#include "FVM/Equation/EquationTerm.hpp"

namespace DREAM {

    class BraamsKarneyRHS : public FVM::EquationTerm {
    public:
        BraamsKarneyRHS(FVM::Grid *grid);

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}

        template<typename F1>
        void SetElements(F1 X);

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        bool SetJacobianBlock(const len_t unknId, const len_t derivId, FVM::Matrix *jac, const real_t *) override;
        void SetMatrixElements(FVM::Matrix *mat, real_t*) override;
        void SetVectorElements(real_t *vec, const real_t *f) override;
    };

}

#endif/*_DREAM_EQUATIONS_BRAAMSKARNEY_RHS_HPP*/
