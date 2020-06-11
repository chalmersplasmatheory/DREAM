#ifndef _DREAM_FVM_DIAGONAL_TERM_HPP
#define _DREAM_FVM_DIAGONAL_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


namespace DREAM::FVM {
    class DiagonalTerm : public EquationTerm {
    private:    
        virtual void DeallocateWeights();

    protected:
        bool hasBeenInitialized = false;
        virtual void AllocateWeights();
        real_t *weights = nullptr;

        virtual void InitializeWeights();
        virtual bool TermDependsOnUnknowns() = 0;
        virtual void SetWeights() = 0;
        virtual void AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) = 0;
    public:
        DiagonalTerm(Grid*);
        ~DiagonalTerm();
        
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        /**
         * This term shows up together with 'PredeterminedParameter' and
         * such, and so we never actually want to assign anything to the
         * vector when evaluating this term (this term indicates that we
         * want to evaluate EVERYTHING ELSE in the equation). */
//        virtual real_t Evaluate(real_t*, const real_t*, const len_t, const len_t);

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual bool GridRebuilt() override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        //virtual void SetMatrixElements(Matrix*, real_t*) override;
        //virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_DIAGONAL_TERM_HPP*/
