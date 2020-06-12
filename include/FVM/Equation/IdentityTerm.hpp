/**
 * Implementation of a term which represents scaleFactor 
 * times the quantity it is applied to.
 */
#ifndef _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP
#define _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP

#include "FVM/Equation/DiagonalLinearTerm.hpp"

namespace DREAM::FVM {
    class IdentityTerm : public DiagonalLinearTerm {
    private:
        real_t scaleFactor;
    protected:
        virtual void SetWeights() override {
            for(len_t i = 0; i<grid->GetNCells(); i++)
                weights[i] = scaleFactor;
        }

    public:
        IdentityTerm(Grid* g, const real_t scaleFactor=1.0) 
            : DiagonalLinearTerm(g), scaleFactor(scaleFactor) {}

    };
}

#endif/*_DREAM_FVM_EQUATION_IDENTITY_TERM_HPP*/
