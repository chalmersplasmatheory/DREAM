/**
 * Implementation of an interpolation object specifically designed
 * to handle the time-evolution of radial profiles of ion charge
 * state densities. The intent of this object is to interpolate the
 * ion charge state density profiles in time.
 */

#include "DREAM/IonInterpolator1D.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
IonInterpolator1D::IonInterpolator1D(
    const len_t nZ0, const len_t nt, const len_t nr,
    const real_t *t, const real_t *densities,
    enum FVM::Interpolator1D::interp_method meth
) {
    this->nZ0 = nZ0;
    this->nt  = nt;
    this->nr  = nr;
    this->interps = new FVM::Interpolator1D*[nZ0];

    this->t         = t;
    this->densities = densities;

    for (len_t i = 0; i < nZ0; i++) {
        this->interps[i] = new FVM::Interpolator1D(
            nt, nr, t, densities+(i*nt*nr),  meth
        );
    }
}

/**
 * Destrutor.
 */
IonInterpolator1D::~IonInterpolator1D() {
    for (len_t i = 0; i < this->nZ0; i++)
        delete this->interps[i];

    delete [] this->interps;

    delete [] this->t;
    delete [] this->densities;
}

/**
 * Evaluate the radial density of the given ion
 * charge state, at the given time.
 *
 * iZ0:  Ion charge state index to evaluate radial density for.
 * t:    Time at which to evaluate the density profile.
 */
const real_t *IonInterpolator1D::Eval(const len_t iZ0, const real_t t) {
    return this->interps[iZ0]->Eval(t);
}

