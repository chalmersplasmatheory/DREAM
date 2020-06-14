/**
 * Implementation of a term which represents an identity term 
 * multiplied by an arbitrary grid-dependent weight function. 
 *      T = w*x
 * where w are grid-dependent weights and x the unknown quantity.
 * Evaluation of weights must be implemented in derived classes. 
 */

#include "FVM/Equation/Operator.hpp"
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
 * Transform the given input vector 'vec' so as to solve for
 * the value of the unknown quantity operated on by this
 * DiagonalLinearTerm. We imagine that this term is part of
 * an equation of the form
 *
 *   w*x + F(y) = 0
 *
 * where F denotes an arbitrary, evaluable function of any
 * other unknown quantities (i.e. y may NOT be x). The input
 * vector is assumed to contain the value 'F', and we
 * subsequently transform the input vector so that the output
 * represents
 *
 *   vec[i] = -F_i(y) / w_i
 *
 * i.e. the corresponding value of the unknown quantity 'x'.
 *
 * vec: Vector, containing the values corresponding to 'F'
 *      to transform.
 */
void DiagonalLinearTerm::EvaluableTransform(
    real_t *vec
) {
    const len_t N = this->grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] = -vec[i] / weights[i];
}

