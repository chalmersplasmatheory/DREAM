/**
 * Implementation of the 'Equation' class, which represents
 * a single physical equation.
 */

#include <algorithm>
#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/TransientTerm.hpp"


using namespace DREAM::FVM;
using namespace std;


/**
 * Constructor.
 */
Equation::Equation(Grid *grid, enum AdvectionDiffusionTerm::advdiff_interpolation intp)
    : grid(grid), advdiff_interpolationMethod(intp) {
    
}

/**
 * Destructor.
 */
Equation::~Equation() {
    // TODO
}

/**
 * Returns the number of non-zero elements inserted into
 * a linear operator matrix by this equation object.
 */
len_t Equation::GetNumberOfNonZerosPerRow() const {
    len_t nnz = 0;

    if (this->tterm != nullptr) nnz = max(nnz, tterm->GetNumberOfNonZerosPerRow());
    if (this->adterm != nullptr) nnz = max(nnz, adterm->GetNumberOfNonZerosPerRow());
    if (this->prescribed != nullptr) nnz = max(nnz, prescribed->GetNumberOfNonZerosPerRow());

    for (auto it = terms.begin(); it != terms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow());

    // Ignore boundary conditions...
    
    return nnz;
}

/**
 * Returns the number of non-zero elements inserted into
 * a jacobian matrix by this equation object.
 */
len_t Equation::GetNumberOfNonZerosPerRow_jac() const {
    len_t nnz = 0;

    if (this->tterm != nullptr) nnz = max(nnz, tterm->GetNumberOfNonZerosPerRow_jac());
    if (this->adterm != nullptr) nnz = max(nnz, adterm->GetNumberOfNonZerosPerRow_jac());
    if (this->prescribed != nullptr) nnz = max(nnz, prescribed->GetNumberOfNonZerosPerRow_jac());

    for (auto it = terms.begin(); it != terms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow());

    // Ignore boundary conditions...
    
    return nnz;
}

/**
 * Build the coefficients of all terms in this equation.
 *
 * t:  Time for which to build the equation terms.
 * dt: Length of time step to take.
 */
void Equation::RebuildTerms(const real_t t, const real_t dt, UnknownQuantityHandler *uqty) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->Rebuild(t, dt, uqty);

    // Advection-diffusion term
    if (adterm != nullptr)
        adterm->Rebuild(t, dt, uqty);

    // TODO Boundary conditions

    // Transient term
    if (tterm != nullptr)
        tterm->Rebuild(t, dt, uqty);
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

