/**
 * Implementation of a term which represents one times the
 * quantity it is applied to.
 */

#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
IdentityTerm::IdentityTerm(Grid *g) : EquationTerm(g) { }

/**
 * Destructor
 */
IdentityTerm::~IdentityTerm() { }


/**
 * Set a block for this term in the given jacobian matrix.
 */
void IdentityTerm::SetJacobianBlock(const len_t derivId, const len_t uqtyId, Matrix *jac) {
    if (derivId == uqtyId) {
        this->SetMatrixElements(jac, nullptr);
    }
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void IdentityTerm::SetMatrixElements(Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, 1);
}

/**
 * Set function vector for this term.
 */
void IdentityTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] = x[i];
}

