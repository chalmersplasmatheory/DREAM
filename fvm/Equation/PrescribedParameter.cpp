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
    : PredeterminedParameter(g), interp_method(interp) {
}
PrescribedParameter::PrescribedParameter(Grid *g, Interpolator1D *i1d)
    : PredeterminedParameter(g), interp(i1d) {
    
    this->time = interp->GetX();
    this->data = interp->GetY();
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
    delete this->interp;
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
        real_t *tm = new real_t[nt];
        real_t *dt = new real_t[nt*N];

        for (len_t i = 0; i < nt; i++)
            tm[i] = t[i];
        for (len_t i = 0; i < N; i++)
            dt[i] = v[i];

        this->time = tm;
        this->data = dt;
    } else {
        this->time = t;
        this->data = v;
    }

    this->interp = new Interpolator1D(nt, N, this->time, this->data, this->interp_method);
}


/**
 * Rebuild (actually, we just need to update the
 * current time, ahead of building the RHS vector).
 *
 * t: Current simulation time.
 */
void PrescribedParameter::Rebuild(const real_t t, const real_t, UnknownQuantityHandler*) {
    if (this->currentTime == t && this->interpolatedData != nullptr)
        return;

    this->currentTime = t;
    this->interpolatedData = this->interp->Eval(t);
}

/**
 * Set the elements in the matrix and on the RHS corresponding
 * to this quantity.
 *
 * mat: Matrix to set elements in (1 is added to the diagonal)
 * rhs: Right-hand-side. Values will be set to the current value of
 *      this parameter.
 */
void PrescribedParameter::SetMatrixElements(Matrix *mat, real_t *rhs) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, 1.0);
    for (len_t i = 0; i < N; i++)
        rhs[i] = -interpolatedData[i];
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution.
 */
void PrescribedParameter::SetVectorElements(real_t *vec, const real_t *x) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] = x[i] - interpolatedData[i];
}

