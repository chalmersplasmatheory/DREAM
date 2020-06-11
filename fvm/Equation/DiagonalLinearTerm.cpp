/**
 * Implementation of a term which represents an identity term 
 * multiplied by an arbitrary grid-dependent weight function. 
 *      T = w*x
 * where w are grid-dependent weights and x the unknown quantity.
 * Evaluation of weights must be implemented in derived classes. 
 */

#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
DiagonalLinearTerm::DiagonalLinearTerm(Grid *g) : DiagonalTerm(g), EvaluableEquationTerm(g) {
    this->grid = g;
    nr = this->DiagonalTerm::nr;
    n1 = this->DiagonalTerm::n1;
    n2 = this->DiagonalTerm::n2;
}


/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void DiagonalLinearTerm::SetMatrixElements(Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, weights[i]);
}

/**
 * Set function vector for this term.
 */
void DiagonalLinearTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        vec[i] += weights[i] * x[i];
}



/**
 * Evaluate this diagonal term. If the term corresponds to
 * the unknown quantity for which the equation applies, the
 * weight is returned so that the caller can appropriately
 * re-scale the result.
 */
real_t *DiagonalLinearTerm::Evaluate(real_t *vec, const real_t *x, const len_t eqnId, const len_t uqtyId) {
    // The diagonal identity term is implied when evaluating equation terms,
    // and so we shouldn't add it here...
    if (eqnId == uqtyId)
        return weights;

    /*len_t N = this->grid->GetNCells();
    const real_t sf = this->scaleFactor;
    for (len_t i = 0; i < N; i++)
        vec[i] -= sf*x[i];*/
    this->SetVectorElements(vec, x);

    return nullptr;
}
