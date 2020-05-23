/**
 * Implementation of a Euler backward transient term.
 *
 * df   f_{n+1} - f_n
 * -- ~ -------------
 * dt        dt
 *
 */

#include <iostream>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
TransientTerm::TransientTerm(Grid *grid, const len_t unknownId)
    : EquationTerm(grid), unknownId(unknownId) { }

/**
 * Rebuild the transient term.
 *
 * dt: Length of next time step to take.
 */
void TransientTerm::Rebuild(const real_t, const real_t dt, UnknownQuantityHandler *uqty) {
    this->dt = dt;
    this->xn = uqty->GetUnknownDataPrevious(this->unknownId);
}

/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 */
void TransientTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *mat
) {
    if (uqtyId == derivId)
        this->SetMatrixElements(mat, nullptr);
}

/**
 * Set the matrix elements corresponding to this
 * transient term.
 *
 * This term assumes that the linearized matrix equation
 * to solve is of the form
 *
 *   df/dt + Mf = -S
 *
 * where M is a linear matrix operator represented by 'mat',
 * and 'S' is a source term stored in 'rhs'.
 *
 * mat: Matrix to set elements of.
 * rhs: Equation RHS.
 */
void TransientTerm::SetMatrixElements(Matrix *mat, real_t *rhs) {
    const len_t N = grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, 1/this->dt, ADD_VALUES);

    if (rhs != nullptr)
        for (len_t i = 0; i < N; i++)
            rhs[i] += this->xn[i] / this->dt;
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution (unused).
 */
void TransientTerm::SetVectorElements(real_t *vec, const real_t *xnp1) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] += (xnp1[i] - xn[i]) / this->dt;
}

