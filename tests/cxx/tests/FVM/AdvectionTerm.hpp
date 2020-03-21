#ifndef _TQSTESTS_FVM_ADVECTION_TERM_HPP
#define _TQSTESTS_FVM_ADVECTION_TERM_HPP

#include "EquationTerm.hpp"


namespace TQSTESTS::FVM {
    class AdvectionTerm : public EquationTerm {
    public:
        AdvectionTerm(const std::string& name) : EquationTerm(name) {}

        virtual bool CheckConservativity(TQS::FVM::RadialGrid*) override;
        virtual bool Run(bool) override;
    };
}

#endif/*_TQSTESTS_FVM_ADVECTION_TERM_HPP*/
