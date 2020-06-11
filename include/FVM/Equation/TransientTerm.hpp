#ifndef _DREAM_FVM_TRANSIENT_TERM_HPP
#define _DREAM_FVM_TRANSIENT_TERM_HPP

#include "FVM/Equation/LinearTransientTerm.hpp"


namespace DREAM::FVM {
    class TransientTerm : public LinearTransientTerm {
    private:
        real_t dt;

        // ID of differentiated quantity
        len_t unknownId;
        // Differentiated quantity at the previous time step
        real_t *xn;
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
        TransientTerm(Grid* g, const len_t unknownId, real_t scaleFactor = 1.0) 
            : LinearTransientTerm(g,unknownId), unknownId(unknownId), scaleFactor(scaleFactor) {}
        
    };
}

#endif/*_DREAM_FVM_TRANSIENT_TERM_HPP*/
