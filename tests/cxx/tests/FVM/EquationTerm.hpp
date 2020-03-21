#ifndef _TQSTESTS_FVM_EQUATION_TERM_HPP
#define _TQSTESTS_FVM_EQUATION_TERM_HPP

#include <string>
#include "FVM/Matrix.hpp"
#include "UnitTest.hpp"

namespace TQSTESTS::FVM {
    class EquationTerm : public UnitTest {
    public:
        EquationTerm(const std::string& s) : UnitTest(s) {}

        bool CheckConservativity();
        virtual bool CheckConservativity(TQS::FVM::RadialGrid*) = 0;
        bool IsConservative(TQS::FVM::Matrix*, TQS::FVM::RadialGrid*, const real_t tol);
    };
}

#endif/*_TQSTESTS_FVM_EQUATION_TERM_HPP*/
