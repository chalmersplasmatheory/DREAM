#ifndef _DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
/**
 * Implementation of a class which represents the Gamma_compton contribution to the n_re equation.
 * Employs the analytical growth rate calculated by RunawayFluid.
 */
namespace DREAM {
    class ComptonRateTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        real_t scaleFactor;
        real_t *gamma;
        real_t *dGamma;
    protected: 
    public:
        ComptonRateTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, 
                RunawayFluid *ref, real_t scaleFactor = 1.0);
        ~ComptonRateTerm();

        void AllocateGamma();
        void DeallocateGamma();
        
        virtual bool GridRebuilt() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;
        virtual void SetWeights() override;

    };
}

#endif /*_DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP*/
