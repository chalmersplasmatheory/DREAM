#ifndef _DREAMTESTS_FVM_ADVECTION_TERM_HPP
#define _DREAMTESTS_FVM_ADVECTION_TERM_HPP

#include "EquationTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAMTESTS::FVM {
    class AdvectionTerm : public EquationTerm {
    public:
        AdvectionTerm(const std::string& name) : EquationTerm(name) {}

        virtual bool CheckConservativity(DREAM::FVM::Grid*) override;
        virtual bool CheckValue(DREAM::FVM::Grid*) override;
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_ADVECTION_TERM_HPP*/
