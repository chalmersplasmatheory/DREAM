/**
 * Implementation of an Euler backward transient term.
 *
 * df   f_{n+1} - f_n
 * -- ~ -------------
 * dt        dt
 *
 */


#ifndef _DREAM_FVM_TRANSIENT_TERM_HPP
#define _DREAM_FVM_TRANSIENT_TERM_HPP

#include "FVM/Equation/LinearTransientTerm.hpp"

namespace DREAM::FVM {
    class TransientTerm : public LinearTransientTerm {
    private:
        real_t scaleFactor;
    protected:
        virtual void SetWeights() override {
            for(len_t i = 0; i < this->grid->GetNCells(); i++)
                weights[i] = scaleFactor;
        }
    public:
        TransientTerm(Grid* g, const len_t unknownId, real_t scaleFactor = 1.0) 
            : LinearTransientTerm(g,unknownId), scaleFactor(scaleFactor) {} 
    };
}

#endif/*_DREAM_FVM_TRANSIENT_TERM_HPP*/
