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

        std::vector<len_t> derivIds;
        std::vector<len_t> derivNMultiples;

        void ResetJacobianColumn();
        void SetPartialWeights(len_t derivId, len_t nMultiples);
        len_t MaxNMultiple() {
            len_t nMultiples = 0;
            for(len_t it=0; it<derivIds.size(); it++)
                if (derivNMultiples[it]>nMultiples)
                    nMultiples = derivNMultiples[it];
            return nMultiples;
        }
        void AllocateDiffWeights() override;
        void DeallocateDiffWeights();
        void ResetDiffWeights();

        void SetElementsInternal(std::function<void(const len_t, const len_t, const real_t)>);

            
    protected:        
        real_t *diffWeights = nullptr;
        UnknownQuantityHandler *unknowns;

        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual void AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetDiffWeights(len_t, len_t) = 0; 
    public:
        DiagonalComplexTerm(Grid*, UnknownQuantityHandler*, Grid *operandGrid=nullptr);
        ~DiagonalComplexTerm();
        
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        // Adds derivId to list of unknown quantities that contributes to Jacobian of this advection term
        void AddUnknownForJacobian(FVM::UnknownQuantityHandler *u, len_t derivId){
            derivIds.push_back(derivId);
            derivNMultiples.push_back(u->GetUnknown(derivId)->NumberOfMultiples());
        }

        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
            { 
                len_t nnz = this->GetNumberOfNonZerosPerRow(); 
                for(len_t i = 0; i<derivIds.size(); i++)
                    nnz += derivNMultiples[i];
                return nnz;
            }

    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_COMPLEX_TERM_HPP*/
