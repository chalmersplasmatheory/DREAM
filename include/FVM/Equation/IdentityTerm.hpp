#ifndef _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP
#define _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
    class IdentityTerm : public EvaluableEquationTerm {
    private:
        real_t scaleFactor;
    public:
        IdentityTerm(Grid*, const real_t scaleFactor=1.0);
        virtual ~IdentityTerm();

        /**
         * This term shows up together with 'PredeterminedParameter' and
         * such, and so we never actually want to assign anything to the
         * vector when evaluating this term (this term indicates that we
         * want to evaluate EVERYTHING ELSE in the equation). */
        virtual void Evaluate(real_t*, const real_t*, const len_t, const len_t);

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override {}

        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_IDENTITY_TERM_HPP*/
