/**
 * Implementation of a prescribed parameter (i.e. a parameter with
 * a given time evolution).
 */

#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
PrescribedParameter::PrescribedParameter(Grid *g, enum Interpolator1D::interp_method interp)
    : EquationTerm(g), interp_method(interp) {
    
}

/**
 * Destructor.
 */
PrescribedParameter::~PrescribedParameter() {
    if (time != nullptr)
        DeallocateData();
}

/**
 * Deallocation routine.
 */
void PrescribedParameter::DeallocateData() {
    delete [] time;
    delete [] data;

    delete [] this->interp;
}

/**
 * Sets the data to interpolate. Note that both the time
 * and data arrays are copied, and can hence be safely free'd
 * outside of this object at any time.
 *
 * nt:   Number of time points.
 * t:    Time array.
 * v:    Parameter value. This array should have N*nt elements, where
 *       N is the number of cells on the grid specified when constructing
 *       this object. The array should be structured logically as an
 *       nt-by-N array, i.e. with in nt blocks with N values each.
 * copy: If 'true', allocates new memory and copies the contents of the
 *       given arrays. Otherwise, just copies the pointers without allocating
 *       any new memory.
 */
void PrescribedParameter::SetData(const len_t nt, real_t *t, real_t *v, bool copy) {
    if (time != nullptr)
        DeallocateData();

    const len_t N = this->grid->GetNCells();
    this->nt = nt;

    if (copy) {
        this->time = new real_t[nt];
        this->data = new real_t[nt*N];

        for (len_t i = 0; i < nt; i++)
            this->time[i] = t[i];
        for (len_t i = 0; i < N; i++)
            this->data[i] = v[i];

        this->interp = new Interpolator1D(nt, N, t, v, this->interp_method);
    } else {
        this->time = t;
        this->data = v;
    }
}


/**
 * Rebuild (actually, we just need to update the
 * current time, ahead of building the RHS vector).
 *
 * t: Current simulation time.
 */
void PrescribedParameter::Rebuild(const real_t t) {
    this->currentTime = t;
}

/**
 * Set the elements on the RHS corresponding to this quantity.
 */
void PrescribedParameter::SetMatrixElements(Matrix*, real_t *rhs) {
    // TODO
}

