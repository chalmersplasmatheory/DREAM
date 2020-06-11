#ifndef _DREAM_FVM_EQUATION_DIAGONAL_QUADRATIC_TERM_HPP
#define _DREAM_FVM_EQUATION_DIAGONAL_QUADRATIC_TERM_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Equation/DiagonalTerm.hpp"


namespace DREAM::FVM {
    class DiagonalQuadraticTerm : public DiagonalTerm {
    private:
        // ID of differentiated quantity
        len_t wUqtyId;
    protected:
        UnknownQuantityHandler *unknowns;
        len_t wUqtyNMultiples;
        virtual len_t GetNumberOfWeightsElements() override
            {return wUqtyNMultiples * grid->GetNCells();}

        virtual bool TermDependsOnUnknowns() override {return false;}
        virtual void AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override;
    public:
        DiagonalQuadraticTerm(Grid*, const len_t, UnknownQuantityHandler*);
        
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_QUADRATIC_TERM_HPP*/
