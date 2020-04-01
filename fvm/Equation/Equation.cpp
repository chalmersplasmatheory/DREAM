/**
 * Implementation of the 'Equation' class, which represents
 * a single physical equation.
 */

#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/TransientTerm.hpp"


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
 * t:  Time for which to build the equation terms.
 * dt: Length of time step to take.
 */
void Equation::RebuildTerms(const real_t t, const real_t dt) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->Rebuild(t);

    // Advection-diffusion term
    if (adterm != nullptr)
        adterm->Rebuild(t);

    // TODO Boundary conditions

    // Transient term
    if (tterm != nullptr)
        tterm->Rebuild(dt);
}

/**
 * Set the matrix elements in the given matrix in order to
 * represent this equation.
 *
 * t:   Time for which to build the matrix.
 * mat: Matrix to set elements of.
 */
void Equation::SetMatrixElements(Matrix *mat, real_t *rhs) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->SetMatrixElements(mat, rhs);

    // Advection-diffusion term?
    if (adterm != nullptr)
        adterm->SetMatrixElements(mat, rhs);

    // Boundary conditions

    // Transient term (must be called last!)
    tterm->SetMatrixElements(mat, rhs);
}

