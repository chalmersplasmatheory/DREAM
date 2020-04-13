#ifndef _DREAMTESTS_FVM_INTERPOLATOR_3D_HPP
#define _DREAMTESTS_FVM_INTERPOLATOR_3D_HPP

#include <functional>
#include "FVM/Interpolator3D.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class Interpolator3D : public UnitTest {
    public:
        struct griddata {
            len_t nx1, nx2, nx3;
            real_t *x1, *x2, *x3;
            real_t *y;
        };
        struct gridlimits { real_t x1min, x1max, x2min, x2max, x3min, x3max; };

        Interpolator3D(const std::string& name) : UnitTest(name) {}

        struct griddata *GenerateData(
            struct gridlimits&,
            const len_t, const len_t, const len_t,
            std::function<real_t(real_t, real_t, real_t)>&
        );
        bool EvalInterpolator3D(
            DREAM::FVM::Interpolator3D*,
            std::function<real_t(real_t, real_t, real_t)>&
        );
        bool TestInterpolation3D(
            enum DREAM::FVM::Interpolator3D::interp_method,
            std::function<real_t(real_t, real_t, real_t)>
        );
        bool TestLinear();
        //bool TestNearest();

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_INTERPOLATOR_3D_HPP*/
