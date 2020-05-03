/**
 * Implementation of an object that interpolates in ADAS
 * rate coefficients.
 */

#include <cmath>
#include <gsl/gsl_interp.h>
#include "DREAM/ADASRateInterpolator.hpp"


using namespace DREAM;

/**
 * Constructor.
 *
 * Z:      Ion atomic charge.
 * nn:     Number of density points.
 * nT:     Number of temperature points.
 * logn:   Logarithm of plasma density (size nn).
 * logT:   Logarithm of plasma temperature (size nT).
 * coeff:  ADAS rate coefficient to interpolate (size Z*nn*nT).
 * interp: 2D interpolation method to use.
 */
ADASRateInterpolator::ADASRateInterpolator(
    const len_t Z, const len_t nn, const len_t nT,
    const real_t *logn, const real_t *logT, const real_t *coeff,
    const gsl_interp2d_type *interp
) : Z(Z), nn(nn), nT(nT), logn(logn), logT(logT), data(coeff) {

    this->splines = new gsl_spline2d*[Z];
    this->nacc = new gsl_interp_accel*[Z];
    this->Tacc = new gsl_interp_accel*[Z];
    const len_t stride = nn*nT;

    for (len_t i = 0; i < Z; i++) {
        this->splines[i] = gsl_spline2d_alloc(interp, nT, nn);
        gsl_spline2d_init(this->splines[i], this->logT, this->logn, this->data + i*stride, nT, nn);

        this->nacc[i] = gsl_interp_accel_alloc();
        this->Tacc[i] = gsl_interp_accel_alloc();
    }
}

/**
 * Destructor.
 */
ADASRateInterpolator::~ADASRateInterpolator() {
    for (len_t i = 0; i < Z; i++) {
        gsl_spline2d_free(this->splines[i]);
        gsl_interp_accel_free(this->nacc[i]);
        gsl_interp_accel_free(this->Tacc[i]);
    }

    delete [] this->Tacc;
    delete [] this->nacc;
    delete [] this->splines;
}

/**
 * Evaluate the spline.
 *
 * Z0: Ion charge state to evaluate coefficient for.
 * n:  Density.
 * T   Temperature.
 */
real_t ADASRateInterpolator::Eval(const len_t Z0, const real_t n, const real_t T) {
    const len_t idx = Z0-1;
    const real_t lT = log10(T);
    const real_t ln = log10(n);

    return gsl_spline2d_eval(
        this->splines[idx], lT, ln, this->Tacc[idx], this->nacc[idx]
    );
}


