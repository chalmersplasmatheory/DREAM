#ifndef _DREAM_ADAS_RATE_INTERPOLATOR_HPP
#define _DREAM_ADAS_RATE_INTERPOLATOR_HPP

#include <cmath>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include "FVM/config.h"

namespace DREAM {
    class ADASRateInterpolator {
    private:
        len_t Z, nn, nT;
        const real_t *logn, *logT, *data;

        // Natural logarithm of 10
        const real_t LN10 = log(10.0);

        /**
         * ADAS rate coefficients are given for n = Z charge
         * states, implying that either of the end charge states
         * (Z0=0 or Z0=Z) yield null rate coefficients. For
         * ionization coefficients, the coefficient should be
         * zero at Z0=Z, while for recombination, it should be
         * zero at Z0=0. Thus, we have to "shift" the recombination
         * coefficient by one (shiftZ0=true).
         */
        bool shiftZ0=false;

        // Interpolation objects
        gsl_interp_accel **nacc, **Tacc;
        //gsl_spline2d **splines;
        gsl_interp2d **interp_c, **interp_l;

    public:
        ADASRateInterpolator(
            const len_t, const len_t, const len_t,
            const real_t*, const real_t*, const real_t*,
            bool shiftZ0,
            const gsl_interp2d_type *interp=gsl_interp2d_bicubic
        );
        virtual ~ADASRateInterpolator();

        real_t Eval(const len_t Z0, const real_t n, const real_t T);

        real_t Eval_deriv_n(const len_t Z0, const real_t n, const real_t T);
        real_t Eval_deriv_T(const len_t Z0, const real_t n, const real_t T);
    };
}

#endif/*_DREAM_ADAS_RATE_INTERPOLATOR_HPP*/
