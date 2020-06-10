#ifndef _DREAM_FVM_EQUATION_WEIGHTED_IDENTITY_TERM_HPP
#define _DREAM_FVM_EQUATION_WEIGHTED_IDENTITY_TERM_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include <algorithm>

namespace DREAM::FVM {
    class WeightedIdentityTerm : public EvaluableEquationTerm {
    protected:
        real_t *weights = nullptr;
        virtual void SetWeights() = 0;
        virtual bool TermDependsOnUnknowns() = 0;
    public:
        WeightedIdentityTerm(Grid*);
        virtual ~WeightedIdentityTerm();

        /**
         * This term shows up together with 'PredeterminedParameter' and
         * such, and so we never actually want to assign anything to the
         * vector when evaluating this term (this term indicates that we
         * want to evaluate EVERYTHING ELSE in the equation). */
        virtual real_t Evaluate(real_t*, const real_t*, const len_t, const len_t);

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override 
            {if(TermDependsOnUnknowns()) SetWeights();}
        virtual bool GridRebuilt() override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_IDENTITY_TERM_HPP*/
