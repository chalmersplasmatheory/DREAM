/**
 * Implementation of a transient term for the ion equations. This term
 * applies to single ion species.
 */

#include "DREAM/Equations/Fluid/IonTransientTerm.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
IonTransientTerm::IonTransientTerm(FVM::Grid *g, IonHandler *ih, const len_t iIon, const len_t unknownId)
    : IonEquationTerm(g, ih, iIon), unknownId(unknownId) { }

/**
 * Rebuild this equation term.
 *
 * t:    Current simulation time.
 * dt:   Time step to take.
 * uqty: List of unknown quantities.
 */
void IonTransientTerm::Rebuild(
    const real_t, const real_t dt, FVM::UnknownQuantityHandler *uqty
) {
    this->dt = dt;
    this->xn = uqty->GetUnknownDataPrevious(this->unknownId);
}

/**
 * Set elements of the jacobian matrix.
 *
 * derivId: ID of unknown quantity with respect to which differentiation
 *          should be carried out.
 * uqtyId:  ID of unknown quantity to differentiate.
 * jac:     Jacobian matrix to build.
 * x:       Current value of the unknown quantity.
 * iIon:    Index of ion to build jacobian for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonTransientTerm::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac,
    const real_t* /*x*/,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if (derivId == uqtyId)
        this->SetCSMatrixElements(jac, nullptr, iIon, Z0, rOffset);
}

/**
 * Set elements of the linear operator matrix.
 *
 * mat:     Linear operator matrix to set elements of.
 * rhs:     Right-hand-side vector.
 * iIon:    Index of ion to build jacobian for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonTransientTerm::SetCSMatrixElements(
    FVM::Matrix *mat, real_t *rhs,
    const len_t /*iIon*/, const len_t /*Z0*/, const len_t rOffset
) {
    const len_t N = grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(rOffset+i, rOffset+i, -1/this->dt);

    if (rhs != nullptr)
        for (len_t i = 0; i < N; i++)
            rhs[rOffset+i] += this->xn[rOffset+i];
}

/**
 * Set elements of the function vector.
 *
 * vec:     Function vector to set elements of.
 * nions:   Ion densities in the current iteration.
 * iIon:    Index of ion to build jacobian for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonTransientTerm::SetCSVectorElements(
    real_t *vec, const real_t *nions,
    const len_t /*iIon*/, const len_t /*Z0*/, const len_t rOffset
) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[rOffset+i] += -(nions[rOffset+i] - xn[rOffset+i]) / this->dt;
}

