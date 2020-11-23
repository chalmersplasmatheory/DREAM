#ifndef _DREAM_FVM_INTERPOLATOR_3D_HPP
#define _DREAM_FVM_INTERPOLATOR_3D_HPP

#include <gsl/gsl_interp.h>
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class Interpolator3D {
    public:
        enum interp_method {
            INTERP_NEAREST,
            INTERP_LINEAR
        };

        enum momentumgrid_type {
            GRID_PXI,
            GRID_PPARPPERP
        };

    private:
        len_t nx1, nx2, nx3;

        const real_t *x1, *x2, *x3;
        const real_t *y;

        enum momentumgrid_type gridtype;
        enum interp_method method;

        bool ownsArrays = true;

        gsl_interp_accel *acc1, *acc2, *acc3;

        len_t _find_x(const real_t, const len_t, const real_t*, gsl_interp_accel*);
        len_t _find_x1(const real_t x) { return _find_x(x, this->nx1, this->x1, this->acc1); }
        len_t _find_x2(const real_t x) { return _find_x(x, this->nx2, this->x2, this->acc2); }
        len_t _find_x3(const real_t x) { return _find_x(x, this->nx3, this->x3, this->acc3); }
        real_t _eval_nearest(const real_t, const real_t, const real_t);
        real_t _eval_linear(const real_t, const real_t, const real_t);

    public:
        Interpolator3D(
            const len_t, const len_t, const len_t,
            const real_t*, const real_t*, const real_t*,
            const real_t*, enum momentumgrid_type,
            enum interp_method meth=INTERP_LINEAR,
            bool ownsArrays=true
        );
        ~Interpolator3D();

        const real_t *Eval(FVM::Grid*, enum momentumgrid_type, enum fluxGridType fgt=FLUXGRIDTYPE_DISTRIBUTION, real_t *out=nullptr);
        const real_t *Eval(
            const len_t, const len_t, const len_t,
            const real_t*, const real_t*, const real_t*,
            enum momentumgrid_type, real_t *out=nullptr
        );

        enum momentumgrid_type GetGridType() const { return this->gridtype; }
    };

    class Interpolator3DException : public FVMException {
    public:
        template<typename ... Args>
        Interpolator3DException(const std::string& msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("Interpolator3D");
        }
    };
}

#endif/*_DREAM_FVM_INTERPOLATOR_3D_HPP*/
