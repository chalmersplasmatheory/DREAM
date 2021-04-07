#ifndef _DREAM_FVM_EQUATION_SCALAR_LINEAR_TERM_HPP
#define _DREAM_FVM_EQUATION_SCALAR_LINEAR_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
    class ScalarLinearTerm : public EquationTerm {
    private:
        virtual void DeallocateWeights();
        virtual void AllocateWeights();
        

    protected:
        len_t nWeights;
        UnknownQuantityHandler *unknowns;
        Grid *targetGrid;
        len_t uqtyId;
        real_t *weights = nullptr;
        virtual void SetWeights() = 0;

        
    public:
        ScalarLinearTerm(Grid*,Grid*, UnknownQuantityHandler*, const len_t);
        
//        virtual void EvaluableTransform(real_t*) override;

        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        /**
         * The following methods are to be inherited from DiagonalTerm
         */
        virtual len_t GetNumberOfNonZerosPerRow() const override 
            {return nWeights;}
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override
            {return GetNumberOfNonZerosPerRow(); }
        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* x) override;
        virtual bool GridRebuilt() override;

    };
}

#endif/*_DREAM_FVM_EQUATION_SCALAR_LINEAR_TERM_HPP*/
