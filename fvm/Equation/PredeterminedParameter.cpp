
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
PredeterminedParameter::PredeterminedParameter(Grid *g)
    : EquationTerm(g) {

    this->currentData = new real_t[g->GetNCells()];
}

/**
 * Destructor.
 */
PredeterminedParameter::~PredeterminedParameter() {
    delete [] this->currentData;
}


/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * jac:     Jacobian matrix block to populate.
 * x:       Value of the unknown quantity.
 *
 * (This term represents a constant, and since the derivative
 * with respect to anything of a constant is zero, we don't need
 * to do anything).
 */
bool PredeterminedParameter::SetJacobianBlock(
    const len_t /*derivId*/, const len_t /*uqtyId*/, Matrix* /*jac*/, const real_t* /*x*/
) { return false; }

/**
 * Set the elements in the matrix and on the RHS corresponding
 * to this quantity.
 *
 * mat: Matrix to set elements in (1 is added to the diagonal)
 * rhs: Right-hand-side. Values will be set to the current value of
 *      this parameter.
 */
void PredeterminedParameter::SetMatrixElements(Matrix* /*mat*/, real_t *rhs) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        rhs[i] = -currentData[i];
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution.
 */
void PredeterminedParameter::SetVectorElements(real_t *vec, const real_t* /*x*/) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] = currentData[i];
}

