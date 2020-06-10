/**
 * Implementation of a term which represents scaleFactor 
 * times the quantity it is applied to.
 */

#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
IdentityTerm::IdentityTerm(Grid *g, const real_t scaleFactor)
    : EvaluableEquationTerm(g), scaleFactor(scaleFactor) { }

/**
 * Destructor
 */
IdentityTerm::~IdentityTerm() { }


/**
 * Evaluate this identity term. If this term corresponds to
 * the unknown quantity for which the equation applies, the
 * scaling factor is returned so that the caller can appropriately
 * re-scale the result.
 */
real_t IdentityTerm::Evaluate(real_t *vec, const real_t *x, const len_t eqnId, const len_t uqtyId) {
    // The diagonal identity term is implied when evaluating equation terms,
    // and so we shouldn't add it here...
    if (eqnId == uqtyId)
        return this->scaleFactor;

    /*len_t N = this->grid->GetNCells();
    const real_t sf = this->scaleFactor;
    for (len_t i = 0; i < N; i++)
        vec[i] -= sf*x[i];*/
    this->SetVectorElements(vec, x);

    return 1;
}

/**
 * Set a block for this term in the given jacobian matrix.
 */
void IdentityTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* /*x*/
) {
    if (derivId == uqtyId) {
        this->SetMatrixElements(jac, nullptr);
    }
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void IdentityTerm::SetMatrixElements(Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();

    const real_t sf = this->scaleFactor;
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, sf);
}

/**
 * Set function vector for this term.
 */
void IdentityTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->grid->GetNCells();

    const real_t sf = this->scaleFactor;
    for (len_t i = 0; i < N; i++)
        vec[i] += sf * x[i];
}

