#ifndef _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP
#define _DREAM_FVM_EQUATION_IDENTITY_TERM_HPP

#include "FVM/Equation/DiagonalLinearTerm.hpp"

namespace DREAM::FVM {
    class IdentityTerm : public DiagonalLinearTerm {
    private:
        real_t scaleFactor;
    protected:
        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                for(len_t i = 0; i < n1[ir]; i++)
                    for(len_t j = 0; j < n2[ir]; j++)
                        weights[offset + n1[ir]*j + i] = scaleFactor;
                offset += n1[ir]*n2[ir];
            }
        }

    public:
        IdentityTerm(Grid* g, const real_t scaleFactor=1.0) 
            : DiagonalLinearTerm(g), scaleFactor(scaleFactor) {}

    };
}

#endif/*_DREAM_FVM_EQUATION_IDENTITY_TERM_HPP*/
