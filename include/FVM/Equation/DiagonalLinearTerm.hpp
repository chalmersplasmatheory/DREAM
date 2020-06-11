#ifndef _DREAM_FVM_EQUATION_DIAGONAL_LINEAR_TERM_HPP
#define _DREAM_FVM_EQUATION_DIAGONAL_LINEAR_TERM_HPP

#include "FVM/Equation/DiagonalTerm.hpp"
#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include <algorithm>

namespace DREAM::FVM {
    class DiagonalLinearTerm : public DiagonalTerm, public EvaluableEquationTerm {
    protected:
//        real_t *weights = nullptr;
//        virtual void SetWeights() = 0;
        virtual bool TermDependsOnUnknowns() override {return false;}
        virtual void AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override {}
        Grid *grid;
        len_t nr;
        len_t *n1, *n2;
        
    public:
        DiagonalLinearTerm(Grid*);
        
        /**
         * This term shows up together with 'PredeterminedParameter' and
         * such, and so we never actually want to assign anything to the
         * vector when evaluating this term (this term indicates that we
         * want to evaluate EVERYTHING ELSE in the equation). */
        virtual real_t* Evaluate(real_t*, const real_t*, const len_t, const len_t) override;

        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        /**
         * The following methods are to be inherited from DiagonalTerm
         */
        virtual len_t GetNumberOfNonZerosPerRow() const override { return this->DiagonalTerm::GetNumberOfNonZerosPerRow(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return this->DiagonalTerm::GetNumberOfNonZerosPerRow_jac(); }
        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override { this->DiagonalTerm::Rebuild(0,0,nullptr); };
        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* x) override
            {return this->DiagonalTerm::SetJacobianBlock(uqtyId,derivId,jac,x);}

    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_LINEAR_TERM_HPP*/
