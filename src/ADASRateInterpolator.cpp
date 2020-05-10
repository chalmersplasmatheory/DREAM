/**
 * Implementation of an object that interpolates in ADAS
 * rate coefficients.
 */

#include <cmath>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
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
    bool shiftZ0, const gsl_interp2d_type *interp
) : Z(Z), nn(nn), nT(nT), logn(logn), logT(logT), data(coeff), shiftZ0(shiftZ0) {

    this->interp = new gsl_interp2d*[Z];
    this->nacc = new gsl_interp_accel*[Z];
    this->Tacc = new gsl_interp_accel*[Z];
    const len_t stride = nn*nT;

    for (len_t i = 0; i < Z; i++) {
        this->interp[i] = gsl_interp2d_alloc(interp, nT, nn);
        gsl_interp2d_init(this->interp[i], this->logT, this->logn, this->data + i*stride, nT, nn);

        this->nacc[i] = gsl_interp_accel_alloc();
        this->Tacc[i] = gsl_interp_accel_alloc();
    }
}

/**
 * Destructor.
 */
ADASRateInterpolator::~ADASRateInterpolator() {
    for (len_t i = 0; i < Z; i++) {
        gsl_interp2d_free(this->interp[i]);
        gsl_interp_accel_free(this->nacc[i]);
        gsl_interp_accel_free(this->Tacc[i]);
    }

    delete [] this->Tacc;
    delete [] this->nacc;
    delete [] this->interp;
}

/**
 * Evaluate the interpolation object.
 *
 * Z0: Ion charge state to evaluate coefficient for.
 * n:  Density.
 * T   Temperature.
 */
real_t ADASRateInterpolator::Eval(const len_t Z0, const real_t n, const real_t T) {
    if ((shiftZ0 && Z0 == 0) || (!shiftZ0 && Z0 == Z))
        return 0;

    const len_t idx = (shiftZ0 ? Z0-1 : Z0);
    const len_t stride = this->nn*this->nT;
    const real_t lT = log10(T);
    const real_t ln = log10(n);

    // coeff = 10^ADASDATA
    return exp(LN10 * gsl_interp2d_eval_extrap(
        // GSL interpolation 2D object
        this->interp[idx],
        // Input data
        this->logT, this->logn, this->data+idx*stride,
        // Point to evaluate data in
        lT, ln,
        // Accelerator objects (for quicker table lookup)
        this->Tacc[idx], this->nacc[idx]
    ));
}


