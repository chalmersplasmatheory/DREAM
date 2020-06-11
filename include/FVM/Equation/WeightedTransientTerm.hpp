#ifndef _DREAM_FVM_WEIGHTED_TRANSIENT_TERM_HPP
#define _DREAM_FVM_WEIGHTED_TRANSIENT_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


namespace DREAM::FVM {
    class WeightedTransientTerm : public EquationTerm {
    private:
        real_t dt;

        // ID of differentiated quantity
        len_t unknownId;
        // Differentiated quantity at the previous time step
        real_t *xn;
        virtual void DeallocateWeights();
        virtual void AllocateWeights();
        virtual void InitializeWeights();
    protected:
        real_t *weights = nullptr;

        virtual bool TermDependsOnUnknowns() = 0;
        virtual void SetWeights() = 0;
    public:
        WeightedTransientTerm(Grid*, const len_t);
        ~WeightedTransientTerm();
        
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual bool GridRebuilt() override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_WEIGHTED_TRANSIENT_TERM_HPP*/
