#ifndef _DREAM_FVM_EQUATION_WEIGHTED_IDENTITY_TERM_HPP
#define _DREAM_FVM_EQUATION_WEIGHTED_IDENTITY_TERM_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include <algorithm>

namespace DREAM::FVM {
    class WeightedIdentityTerm : public EvaluableEquationTerm {
    private:
        real_t *weights;

        std::function<real_t(len_t,len_t,len_t)> *weightFunc;
    public:
        WeightedIdentityTerm(Grid*,  std::function<real_t(len_t,len_t,len_t)> *weightFun);
        virtual ~WeightedIdentityTerm();

        /**
         * This term shows up together with 'PredeterminedParameter' and
         * such, and so we never actually want to assign anything to the
         * vector when evaluating this term (this term indicates that we
         * want to evaluate EVERYTHING ELSE in the equation). */
        virtual void Evaluate(real_t*, const real_t*, const len_t, const len_t);

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_IDENTITY_TERM_HPP*/
