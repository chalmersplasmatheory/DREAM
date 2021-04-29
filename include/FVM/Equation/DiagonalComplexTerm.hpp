#ifndef _DREAM_FVM_EQUATION_DIAGONAL_COMPLEX_TERM_HPP
#define _DREAM_FVM_EQUATION_DIAGONAL_COMPLEX_TERM_HPP

#include <functional>
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Equation/DiagonalTerm.hpp"


namespace DREAM::FVM {
    class DiagonalComplexTerm : public DiagonalTerm {
    private:
        // Grid on which the operand (i.e. quantity on which this term operates)
        // of this term lives. May be a nullptr, in which case it is assumed to
        // be the same as 'grid'.
        Grid *operandGrid;

        void ResetJacobianColumn();
        void SetPartialWeights(len_t derivId, len_t nMultiples);
        void AllocateDiffWeights() override;
        void DeallocateDiffWeights();
        void ResetDiffWeights();

        void SetElementsInternal(std::function<void(const len_t, const len_t, const real_t)>);

            
    protected:        
        real_t *diffWeights = nullptr;
        UnknownQuantityHandler *unknowns;

        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual bool AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetDiffWeights(len_t, len_t) = 0; 
    public:
        DiagonalComplexTerm(Grid*, UnknownQuantityHandler*, Grid *operandGrid=nullptr);
        ~DiagonalComplexTerm();
        
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_COMPLEX_TERM_HPP*/
