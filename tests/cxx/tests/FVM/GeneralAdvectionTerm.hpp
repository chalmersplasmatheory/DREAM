#ifndef _DREAMTESTS_FVM_EQUATION_GENERAL_ADVECTION_TERM_HPP
#define _DREAMTESTS_FVM_EQUATION_GENERAL_ADVECTION_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "AdvectionTerm.hpp"

namespace DREAMTESTS::FVM {
    class GeneralAdvectionTerm
        : public DREAM::FVM::AdvectionTerm {
    public:
        GeneralAdvectionTerm(DREAM::FVM::Grid*);
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAMTESTS_FVM_EQUATION_GENERAL_ADVECTION_TERM_HPP*/
