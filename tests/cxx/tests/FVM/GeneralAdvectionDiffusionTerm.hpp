#ifndef _DREAMTESTS_FVM_EQUATION_GENERAL_ADVECTION_DIFFUSION_TERM_HPP
#define _DREAMTESTS_FVM_EQUATION_GENERAL_ADVECTION_DIFFUSION_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "AdvectionTerm.hpp"

namespace DREAMTESTS::FVM {
    class GeneralAdvectionDiffusionTerm
        : public DREAM::FVM::AdvectionDiffusionTerm {
    public:
        GeneralAdvectionDiffusionTerm(DREAM::FVM::Grid*);
    };
}

#endif/*_DREAMTESTS_FVM_EQUATION_GENERAL_ADVECTION_DIFFUSION_TERM_HPP*/
