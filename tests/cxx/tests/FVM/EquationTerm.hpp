#ifndef _DREAMTESTS_FVM_EQUATION_TERM_HPP
#define _DREAMTESTS_FVM_EQUATION_TERM_HPP

#include <string>
#include "FVM/Matrix.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class EquationTerm : public UnitTest {
    public:
        EquationTerm(const std::string& s) : UnitTest(s) {}

        bool CheckConservativity();
        virtual bool CheckConservativity(DREAM::FVM::RadialGrid*) = 0;
        bool IsConservative(DREAM::FVM::Matrix*, DREAM::FVM::RadialGrid*, const real_t tol);
    };
}

#endif/*_DREAMTESTS_FVM_EQUATION_TERM_HPP*/
