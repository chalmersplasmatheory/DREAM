/**
 * Implementation of the 'Equation' class, which represents
 * a single physical equation.
 */

#include <algorithm>
#include "FVM/Equation/Equation.hpp"
//#include "FVM/Equation/TransientTerm.hpp"


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

    //if (this->tterm != nullptr) nnz = max(nnz, tterm->GetNumberOfNonZerosPerRow());
    if (this->adterm != nullptr) nnz = max(nnz, adterm->GetNumberOfNonZerosPerRow());
    if (this->predetermined != nullptr) nnz = max(nnz, predetermined->GetNumberOfNonZerosPerRow());

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

    //if (this->tterm != nullptr) nnz = max(nnz, tterm->GetNumberOfNonZerosPerRow_jac());
    if (this->adterm != nullptr) nnz = max(nnz, adterm->GetNumberOfNonZerosPerRow_jac());
    if (this->predetermined != nullptr) nnz = max(nnz, predetermined->GetNumberOfNonZerosPerRow_jac());

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
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->Rebuild(t, uqty);

    // Transient term
    /*if (tterm != nullptr)
        tterm->Rebuild(t, dt, uqty);*/
}

/**
 * Set the specified block in the given jacobian matrix.
 *
 * uqtyId:  ID of the unknown quantity to which the matrix row belongs.
 * derivId: ID of the unknown quantity with respect to which the
 *          equation should be differentiated.
 * jac:     Jacobian matrix (block) to set.
 */
void Equation::SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->SetJacobianBlock(derivId, uqtyId, jac);

    // Advection-diffusion term?
    if (adterm != nullptr)
        adterm->SetJacobianBlock(derivId, uqtyId, jac);

    // TODO Boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->AddToJacobianBlock(derivId, uqtyId, jac);

}

/**
 * Routine specifically designed for boundary conditions which need to
 * overwrite elements in the jacobian matrix.
 *
 * uqtyId:  ID of the unknown quantity to which the matrix row belongs.
 * derivId: ID of the unknown quantity with respect to which the
 *          equation should be differentiated.
 * jac:     Jacobian matrix (block) to set.
 */
void Equation::SetJacobianBlockBC(const len_t uqtyId, const len_t derivId, Matrix *jac) {
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        adterm->SetJacobianBlock(derivId, uqtyId, jac);
}

/**
 * Set the linear operator matrix elements in the given
 * matrix in order to represent this equation.
 *
 * mat: Matrix to set elements of.
 * rhs: Vector representing equation right-hand-side.
 */
void Equation::SetMatrixElements(Matrix *mat, real_t *rhs) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->SetMatrixElements(mat, rhs);

    // Advection-diffusion term?
    if (adterm != nullptr)
        adterm->SetMatrixElements(mat, rhs);

    // Boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->AddToMatrixElements(mat, rhs);

    // TODO Partially assemble matrix
    mat->PartialAssemble();

    // Hard set boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->SetMatrixElements(mat, rhs);
}

/**
 * Evaluate this equation and assign its value to the
 * given function vector.
 * 
 * vec: Function vector to assign evaluated equation to.
 * x:   Value of the unknown to evaluate the function for.
 */
void Equation::SetVectorElements(real_t *vec, const real_t *x) {
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->SetVectorElements(vec, x);

    // Advection-diffusion term?
    if (adterm != nullptr)
        adterm->SetVectorElements(vec, x);

    // Boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->AddToVectorElements(vec, x);

    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->SetVectorElements(vec, x);
}

