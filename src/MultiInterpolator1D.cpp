/**
 * Implementation of an interpolation object specifically designed
 * to handle the time-evolution of radial profiles of ion charge
 * state densities. The intent of this object is to interpolate the
 * ion charge state density profiles in time.
 */

#include "DREAM/MultiInterpolator1D.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
MultiInterpolator1D::MultiInterpolator1D(
    const len_t nZ0, const len_t nt, const len_t nr,
    const real_t *t, const real_t *x,
    enum FVM::Interpolator1D::interp_method meth
) {
    this->nZ0 = nZ0;
    this->nt  = nt;
    this->nr  = nr;
    this->interps = new FVM::Interpolator1D*[nZ0];

    this->t = t;
    this->x = x;

    for (len_t i = 0; i < nZ0; i++) {
        this->interps[i] = new FVM::Interpolator1D(
            nt, nr, t, x+(i*nt*nr),  meth,
			// owns the data "x+(i*nt*nr)"? (i.e. should it be deleted
			// by this Interpolator1D object?)
			false
        );
    }
}

/**
 * Destrutor.
 */
MultiInterpolator1D::~MultiInterpolator1D() {
    for (len_t i = 0; i < this->nZ0; i++)
        delete this->interps[i];

    delete [] this->interps;

    delete [] this->t;
    delete [] this->x;
}

/**
 * Evaluate the radial density of the given parameter
 * index, at the given time.
 *
 * n: Index to evaluate radial density for.
 * t: Time at which to evaluate the density profile.
 */
const real_t *MultiInterpolator1D::Eval(const len_t n, const real_t t) {
    return this->interps[n]->Eval(t);
}

