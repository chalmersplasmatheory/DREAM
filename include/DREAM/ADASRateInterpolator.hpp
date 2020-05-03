#ifndef _DREAM_ADAS_RATE_INTERPOLATOR_HPP
#define _DREAM_ADAS_RATE_INTERPOLATOR_HPP

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline2d.h>
#include "FVM/config.h"

namespace DREAM {
    class ADASRateInterpolator {
    private:
        len_t Z, nn, nT;
        const real_t *logn, *logT, *data;

        // Interpolation objects
        gsl_interp_accel **nacc, **Tacc;
        gsl_spline2d **splines;

    public:
        ADASRateInterpolator(
            const len_t, const len_t, const len_t,
            const real_t*, const real_t*, const real_t*,
            const gsl_interp2d_type *interp=gsl_interp2d_bicubic
        );
        virtual ~ADASRateInterpolator();

        real_t Eval(const len_t, const real_t, const real_t);
    };
}

#endif/*_DREAM_ADAS_RATE_INTERPOLATOR_HPP*/
