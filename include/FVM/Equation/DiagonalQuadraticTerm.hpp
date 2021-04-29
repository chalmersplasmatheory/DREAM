#ifndef _DREAM_FVM_EQUATION_DIAGONAL_QUADRATIC_TERM_HPP
#define _DREAM_FVM_EQUATION_DIAGONAL_QUADRATIC_TERM_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Equation/DiagonalTerm.hpp"


namespace DREAM::FVM {
    class DiagonalQuadraticTerm : public DiagonalTerm, public EvaluableEquationTerm {
    private:
        // ID of differentiated quantity
        len_t wUqtyId;
    protected:
        FVM::Grid *grid;
        len_t nr, *n1, *n2;

        UnknownQuantityHandler *unknowns;
        len_t wUqtyNMultiples;
        virtual len_t GetNumberOfWeightsElements() override
            {return wUqtyNMultiples * this->grid->GetNCells();}

        virtual bool TermDependsOnUnknowns() override {return false;}
        virtual bool AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override;
    public:
        DiagonalQuadraticTerm(Grid*, const len_t, UnknownQuantityHandler*);
        
        virtual len_t GetNumberOfNonZerosPerRow() const override { return this->DiagonalTerm::GetNumberOfNonZerosPerRow(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
            { return this->DiagonalTerm::GetNumberOfNonZerosPerRow_jac() + this->wUqtyNMultiples; }

        virtual void EvaluableTransform(real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override { this->DiagonalTerm::Rebuild(0,0,nullptr); };
        virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* x) override
            {return this->DiagonalTerm::SetJacobianBlock(uqtyId,derivId,jac,x);}
    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_QUADRATIC_TERM_HPP*/
