#ifndef _DREAMTESTS_FVM_EQUATION_TERM_HPP
#define _DREAMTESTS_FVM_EQUATION_TERM_HPP

#include <functional>
#include <string>
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class EquationTerm : public UnitTest {
    public:
        EquationTerm(const std::string& s) : UnitTest(s) {}

        bool CheckConservativity();
        virtual bool CheckConservativity(DREAM::FVM::Grid*) = 0;
        bool CheckValue();
        virtual bool CheckValue(DREAM::FVM::Grid*) = 0;

        void EvaluateFVM(
            const len_t, const len_t, const len_t,
            const real_t*, const real_t*, const real_t*,
            std::function<real_t(len_t, len_t, len_t)>,    // Vp
            std::function<real_t(len_t, len_t, len_t)>,    // Vp_fr
            std::function<real_t(len_t, len_t, len_t)>,    // Vp_f1
            std::function<real_t(len_t, len_t, len_t)>,    // Vp_f2
            std::function<real_t(len_t, len_t, len_t)>,    // fluxR
            std::function<real_t(len_t, len_t, len_t)>,    // flux1
            std::function<real_t(len_t, len_t, len_t)>,    // flux2
            real_t*
        );
        bool IsConservative(DREAM::FVM::Matrix*, DREAM::FVM::Grid*, const real_t tol);
        bool IsReallyConservative(DREAM::FVM::Matrix*, DREAM::FVM::Grid*, const real_t tol);
    };
}

#endif/*_DREAMTESTS_FVM_EQUATION_TERM_HPP*/
