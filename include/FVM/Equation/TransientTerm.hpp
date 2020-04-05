#ifndef _DREAM_FVM_TRANSIENT_TERM_HPP
#define _DREAM_FVM_TRANSIENT_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


namespace DREAM::FVM {
    class TransientTerm : public EquationTerm {
    private:
        real_t dt;

        // ID of differentiated quantity
        len_t unknownId;
        // Differentiated quantity at the previous time step
        real_t *xn;

    public:
        TransientTerm(Grid*, const len_t);
        
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_TRANSIENT_TERM_HPP*/
