#ifndef _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP
#define _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
    class IdentityTerm : public EquationTerm {
    public:
        IdentityTerm(Grid*);
        virtual ~IdentityTerm();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override {}

        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_IDENTITY_TERM_HPP*/
