#ifndef _DREAMTESTS_FVM_INTERPOLATOR_1D_HPP
#define _DREAMTESTS_FVM_INTERPOLATOR_1D_HPP

#include <functional>
#include "FVM/Interpolator1D.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class Interpolator1D : public UnitTest {
    public:
        Interpolator1D(const std::string& name) : UnitTest(name) {}

        bool TestLinear_general(
            enum DREAM::FVM::Interpolator1D::interp_method,
            std::function<bool(len_t, len_t, real_t, const real_t*, const real_t*, const real_t*, const std::function<real_t(const real_t, const real_t)>&)>);
        bool TestLinear();
        bool TestNearest();

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_INTERPOLATOR_1D_HPP*/
