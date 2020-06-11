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
        UnknownQuantityHandler *unknowns;
        // ID of differentiated quantity
        len_t wUqtyId;

    protected:
        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual void AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override;
    public:
        DiagonalQuadraticTerm(Grid*, const len_t, UnknownQuantityHandler*);
        ~DiagonalQuadraticTerm();
        
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_QUADRATIC_TERM_HPP*/
