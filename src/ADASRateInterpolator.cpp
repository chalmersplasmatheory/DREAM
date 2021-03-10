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

    this->interp_c = new gsl_interp2d*[Z];
    this->interp_l = new gsl_interp2d*[Z];
    this->nacc = new gsl_interp_accel*[Z];
    this->Tacc = new gsl_interp_accel*[Z];
    const len_t stride = nn*nT;

    for (len_t i = 0; i < Z; i++) {
        this->interp_c[i] = gsl_interp2d_alloc(interp, nn, nT);
        gsl_interp2d_init(this->interp_c[i], this->logn, this->logT, this->data + i*stride, nn, nT);

        if (interp == gsl_interp2d_bilinear)
            this->interp_l[i] = this->interp_c[i];
        else {
            this->interp_l[i] = gsl_interp2d_alloc(gsl_interp2d_bilinear, nn, nT);
            gsl_interp2d_init(this->interp_l[i], this->logn, this->logT, this->data + i*stride, nn, nT);
        }

        this->nacc[i] = gsl_interp_accel_alloc();
        this->Tacc[i] = gsl_interp_accel_alloc();
    }
}

/**
 * Destructor.
 */
ADASRateInterpolator::~ADASRateInterpolator() {
    for (len_t i = 0; i < Z; i++) {
        if (this->interp_c[i] != this->interp_l[i])
            gsl_interp2d_free(this->interp_l[i]);

        gsl_interp2d_free(this->interp_c[i]);
        gsl_interp_accel_free(this->nacc[i]);
        gsl_interp_accel_free(this->Tacc[i]);
    }

    delete [] this->Tacc;
    delete [] this->nacc;
    delete [] this->interp_l;
    delete [] this->interp_c;
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

    // presumably these are expected to be 0 within numerical errors.
    if((n<=0) || (T<=0))
        return 0;

    const len_t idx = (shiftZ0 ? Z0-1 : Z0);
    const len_t stride = this->nn*this->nT;
    const real_t lT = log10(T);
    const real_t ln = log10(n);

    // coeff = 10^ADASDATA
    gsl_interp2d *spln;
    if (ln < this->logn[0] || lT < this->logT[0] || ln > this->logn[this->nn-1] || lT > this->logT[this->nT-1])
        spln = this->interp_l[idx];
    else
        spln = this->interp_c[idx];

    return exp(LN10 * gsl_interp2d_eval_extrap(
        // GSL interpolation 2D object
        spln,
        // Input data
        this->logn, this->logT, this->data+idx*stride,
        // Point to evaluate data in
        ln, lT,
        // Accelerator objects (for quicker table lookup)
        this->nacc[idx], this->Tacc[idx]
    ));
}

real_t ADASRateInterpolator::Eval_deriv_n(const len_t Z0, const real_t n, const real_t T) {
    real_t eps = std::numeric_limits<real_t>::epsilon();
    real_t dn = sqrt(eps)*n + eps;
    return ( Eval(Z0,n+dn,T) - Eval(Z0,n,T) ) / dn;
}


real_t ADASRateInterpolator::Eval_deriv_T(const len_t Z0, const real_t n, const real_t T) {
    real_t eps = std::numeric_limits<real_t>::epsilon();
    real_t dT = sqrt(eps)*T + eps;
    return ( Eval(Z0,n,T+dT) - Eval(Z0,n,T) ) / dT;
}



