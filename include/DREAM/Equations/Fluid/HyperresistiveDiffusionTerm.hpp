#ifndef _DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class HyperresistiveDiffusionTerm
        : public FVM::DiffusionTerm {
    private:
    real_t *Lambda; 
    real_t *psi_t;
    public:
        HyperresistiveDiffusionTerm(FVM::Grid*, real_t*, real_t*);
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_POLOIDAL_FLUX_HYPERRESISTIVE_DIFFUSION_TERM_HPP*/


