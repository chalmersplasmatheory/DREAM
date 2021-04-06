#ifndef _DREAM_FVM_EQUATION_DIAGONAL_LINEAR_TERM_HPP
#define _DREAM_FVM_EQUATION_DIAGONAL_LINEAR_TERM_HPP

#include "FVM/Equation/DiagonalTerm.hpp"
#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include <algorithm>

namespace DREAM::FVM {
    class DiagonalLinearTerm : public DiagonalTerm, public EvaluableEquationTerm {
    protected:
        virtual bool TermDependsOnUnknowns() override {return false;}
        virtual bool AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override {return false;}
        Grid *grid;
        len_t nr;
        len_t *n1, *n2;

        len_t nMultiples;
        
        
    public:
        DiagonalLinearTerm(Grid*);
        
        virtual void EvaluableTransform(real_t*) override;

        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        /**
         * The following methods are to be inherited from DiagonalTerm
         */
        virtual len_t GetNumberOfNonZerosPerRow() const override { return this->DiagonalTerm::GetNumberOfNonZerosPerRow(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return this->DiagonalTerm::GetNumberOfNonZerosPerRow_jac(); }
        virtual void Rebuild(const real_t t, const real_t dt, UnknownQuantityHandler *u) override { this->DiagonalTerm::Rebuild(t,dt,u); };
        virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* x) override
            {return this->DiagonalTerm::SetJacobianBlock(uqtyId,derivId,jac,x);}

    };
}

#endif/*_DREAM_FVM_EQUATION_DIAGONAL_LINEAR_TERM_HPP*/
