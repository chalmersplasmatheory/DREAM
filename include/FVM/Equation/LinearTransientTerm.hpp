#ifndef _DREAM_FVM_EQUATION_LINEAR_TRANSIENT_TERM_HPP
#define _DREAM_FVM_EQUATION_LINEAR_TRANSIENT_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/DiagonalTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


namespace DREAM::FVM {
    class LinearTransientTerm : public DiagonalTerm {
    private:
        real_t dt;

        // ID of differentiated quantity
        len_t unknownId;
        // Differentiated quantity at the previous time step
        real_t *xn;
    protected:
        virtual bool TermDependsOnUnknowns() override {return false;}
        virtual bool AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) override {return false;}
    public:
        LinearTransientTerm(Grid*, const len_t);
        
        virtual void Rebuild(const real_t, const real_t dt, UnknownQuantityHandler *uqty) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_LINEAR_TRANSIENT_TERM_HPP*/
