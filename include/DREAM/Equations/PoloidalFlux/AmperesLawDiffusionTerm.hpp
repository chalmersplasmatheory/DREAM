#ifndef _DREAM_EQUATIONS_POLOIDAL_FLUX_AMPERES_LAW_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_POLOIDAL_FLUX_AMPERES_LAW_DIFFUSION_TERM_HPP

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class AmperesLawDiffusionTerm
        : public FVM::DiffusionTerm {
    private:
        //virtual void SetPartialDiffusionTerm(len_t, len_t) override{};
    public:
        AmperesLawDiffusionTerm(FVM::Grid*);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_POLOIDAL_FLUX_AMPERES_LAW_DIFFUSION_TERM_HPP*/


