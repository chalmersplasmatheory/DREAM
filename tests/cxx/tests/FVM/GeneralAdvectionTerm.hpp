#ifndef _TQSTESTS_FVM_EQUATION_GENERAL_ADVECTION_TERM_HPP
#define _TQSTESTS_FVM_EQUATION_GENERAL_ADVECTION_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "AdvectionTerm.hpp"

namespace TQSTESTS::FVM {
    class GeneralAdvectionTerm
        : public TQS::FVM::AdvectionTerm {
    public:
        GeneralAdvectionTerm(TQS::FVM::RadialGrid*);
        virtual void Rebuild(const real_t) override;
    };
}

#endif/*_TQSTESTS_FVM_EQUATION_GENERAL_ADVECTION_TERM_HPP*/
