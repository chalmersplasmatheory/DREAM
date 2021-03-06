#ifndef _DREAMTESTS_FVM_EQUATION_GENERAL_DIFFUSION_TERM_HPP
#define _DREAMTESTS_FVM_EQUATION_GENERAL_DIFFUSION_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DiffusionTerm.hpp"

namespace DREAMTESTS::FVM {
    class GeneralDiffusionTerm
        : public DREAM::FVM::DiffusionTerm {
    private:
        real_t value;

    public:
        GeneralDiffusionTerm(DREAM::FVM::Grid*, const real_t v=0.0);
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAMTESTS_FVM_EQUATION_GENERAL_DIFFUSION_TERM_HPP*/

