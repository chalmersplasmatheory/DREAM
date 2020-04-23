#ifndef _DREAM_ION_INTERPOLATOR_1D_HPP
#define _DREAM_ION_INTERPOLATOR_1D_HPP

#include "FVM/Interpolator1D.hpp"

namespace DREAM {
    class IonInterpolator1D {
    private:
        len_t nZ0=0, nt, nr;
        FVM::Interpolator1D **interps=nullptr;

        const real_t *t, *densities;

    public:
        IonInterpolator1D(
            const len_t, const len_t, const len_t,
            const real_t*, const real_t*,
            enum FVM::Interpolator1D::interp_method
        );
        ~IonInterpolator1D();

        const real_t *Eval(const len_t, const real_t);
    };
}

#endif/*_DREAM_ION_INTERPOLATOR_1D_HPP*/
