#ifndef _DREAM_MULTI_INTERPOLATOR_1D_HPP
#define _DREAM_MULTI_INTERPOLATOR_1D_HPP

#include "FVM/Interpolator1D.hpp"

namespace DREAM {
    class MultiInterpolator1D {
    private:
        len_t nZ0=0, nt, nr;
        FVM::Interpolator1D **interps=nullptr;

        const real_t *t, *x;

    public:
        MultiInterpolator1D(
            const len_t, const len_t, const len_t,
            const real_t*, const real_t*,
            enum FVM::Interpolator1D::interp_method
        );
        ~MultiInterpolator1D();

        const real_t *Eval(const len_t, const real_t);
    };
}

#endif/*_DREAM_MULTI_INTERPOLATOR_1D_HPP*/
