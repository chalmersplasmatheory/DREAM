#ifndef _DREAMTESTS_FVM_DIFFUSION_TERM_HPP
#define _DREAMTESTS_FVM_DIFFUSION_TERM_HPP

#include "EquationTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAMTESTS::FVM {
    class DiffusionTerm : public EquationTerm {
    public:
        DiffusionTerm(const std::string& name) : EquationTerm(name) {}

        virtual bool CheckConservativity(DREAM::FVM::Grid*) override;
        virtual bool CheckValue(DREAM::FVM::Grid*) override;
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_DIFFUSION_TERM_HPP*/
