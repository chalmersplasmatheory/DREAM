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

    public:
        TransientTerm(Grid*);
        
        virtual void Rebuild(const real_t) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_TRANSIENT_TERM_HPP*/
