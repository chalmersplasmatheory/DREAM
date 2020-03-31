/**
 * Implementation of the 'Equation' class, which represents
 * a single physical equation.
 */

#include "FVM/Equation/Equation.hpp"


using namespace DREAM::FVM;


/**
 * Destructor.
 */
Equation::~Equation() {
    // TODO
}

/**
 * Build the coefficients of all terms in this equation.
 *
 * t: Time for which to build the equation terms.
 */
void Equation::RebuildTerms(const real_t t) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->Rebuild(t);
}

/**
 * Set the matrix elements in the given matrix in order to
 * represent this equation.
 *
 * t:   Time for which to build the matrix.
 * mat: Matrix to set elements of.
 */
void Equation::SetMatrixElements(Matrix *mat) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->SetMatrixElements(mat);
}

