#ifndef _DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class ComptonRateTerm : public FVM::DiagonalComplexTerm, public RunawaySourceTerm {
    private:
        RunawayFluid *REFluid;
        real_t scaleFactor;
        real_t *dGamma = nullptr;
        len_t nr_tmp = 0;
        void AllocateDGamma();
    public:
        ComptonRateTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, 
                RunawayFluid *ref, FVM::Grid*, real_t scaleFactor = 1.0);
        ~ComptonRateTerm();
        
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;
        virtual void SetWeights() override;

    };
}

#endif /*_DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP*/
