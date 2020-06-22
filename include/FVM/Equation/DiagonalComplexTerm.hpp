#ifndef _DREAM_FVM_EQUATION_DIAGONAL_COMPLEX_TERM_HPP
#define _DREAM_FVM_EQUATION_DIAGONAL_COMPLEX_TERM_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Equation/DiagonalTerm.hpp"


namespace DREAM::FVM {
    class DiagonalComplexTerm : public DiagonalTerm {
    private:

        std::vector<len_t> derivIds;
        std::vector<len_t> derivNMultiples;

        void ResetJacobianColumn();
        void SetPartialWeights(len_t derivId, len_t nMultiples);
        len_t MaxNMultiple()
            {
            len_t nMultiples = 0;
            for(len_t it=0; it<derivIds.size(); it++)
                if (derivNMultiples[it]>nMultiples)
                    nMultiples = derivNMultiples[it];
            return nMultiples;
            }
            void AllocateDiffWeights();
            void DeallocateDiffWeights();
            void ResetDiffWeights();
            
    protected:        
        real_t *diffWeights = nullptr;
        UnknownQuantityHandler *unknowns;

        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual void AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetDiffWeights(len_t, len_t) = 0; 
    public:
        DiagonalComplexTerm(Grid*, UnknownQuantityHandler*);
        ~DiagonalComplexTerm();
        
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        // Adds derivId to list of unknown quantities that contributes to Jacobian of this advection term
        void AddUnknownForJacobian(FVM::UnknownQuantityHandler *u, len_t derivId){
            derivIds.push_back(derivId);
            derivNMultiples.push_back(u->GetUnknown(derivId)->NumberOfMultiples());
        }

    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_COMPLEX_TERM_HPP*/
