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
 * Evaluate the terms of this equation.
 *
 * vec: Vector to store value in.
 * x:   Current value of the unknown quantity to which this
 *      operator is applied.
 */
void Equation::Evaluate(real_t *vec, const real_t *x) {
    if (IsPredetermined()) {
        const real_t *data = this->predetermined->GetData();
        for (len_t i = 0; i < this->grid->GetNCells(); i++)
            vec[i] += data[i];
    } else {
        this->SetVectorElements(vec, x);
    }
}

/**
 * If this equation contains a single EvaluableTerm, applies
 * the transform of that EvaluableTerm to the given vector.
 *
 * Essentially, it is assumed that this equation represents
 * a single (invertible) operator f(x), operating on the
 * unknown quantity 'x'. Further, it is assumed that this term
 * is part of an equation of the form
 *
 *   f(x) + g(y) = 0
 *
 * where 'y' represents one or more unknown quantities, NOT
 * INCLUDING x. Given the vector 'vec', which is assumed to
 * contain 'g(y)', this method calculates
 *
 *   x = f^-1( -g(y) ),
 *
 * where 'f^-1' denotes the inverse of 'f(x)'.
 */
void Equation::EvaluableTransform(real_t *vec) {
    if (this->eval_terms.size() != 1)
        throw EquationException(
            "This equation must have exactly one evaluable term for it to be evaluable."
        );

    eval_terms[0]->EvaluableTransform(vec);
}

/**
 * Returns true if this equation is evaluable, i.e. if it
 * consists of exactly one EvaluableEquationTerm. This means
 * that we can solve for the unknown quantity to which this
 * operator is applied in the manner described in the comment
 * for 'EvaluableTransform()' above.
 *
 * We also return 'true' for 'PredeterminedParameter's, since
 * these have an implied IdentityTerm, which is evaluable, in them.
 */
bool Equation::IsEvaluable() const {
    if (this->adterm != nullptr ||
        boundaryConditions.size() > 0 ||
        terms.size() > 0)
        return false;
    else return
        (this->predetermined != nullptr && this->eval_terms.size() == 0) ||
        (this->predetermined == nullptr && this->eval_terms.size() == 1);
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
    for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
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
    for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
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
    // Predetermined value
    if (predetermined != nullptr) {
        this->predetermined->Rebuild(t, dt, uqty);
        return;
    }

    // Evaluatable equation terms
    for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
        (*it)->Rebuild(t, dt, uqty);

    // Other equation terms
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->Rebuild(t, dt, uqty);

    // Advection-diffusion term
    if (adterm != nullptr)
        adterm->Rebuild(t, dt, uqty);

    // Boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->Rebuild(t, uqty);
}

/**
 * Set the specified block in the given jacobian matrix.
 *
 * uqtyId:  ID of the unknown quantity to which the matrix row belongs.
 * derivId: ID of the unknown quantity with respect to which the
 *          equation should be differentiated.
 * jac:     Jacobian matrix (block) to set.
 * x:       Value of the unknown quantity.
 */
void Equation::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x
) {
    for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
        (*it)->SetJacobianBlock( uqtyId, derivId, jac, x);

    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);

    // Advection-diffusion term?
    if (adterm != nullptr)
        adterm->SetJacobianBlock(uqtyId, derivId, jac, x);

    // Boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->AddToJacobianBlock(uqtyId, derivId, jac, x);
}

/**
 * Routine specifically designed for boundary conditions which need to
 * overwrite elements in the jacobian matrix.
 *
 * uqtyId:  ID of the unknown quantity to which the matrix row belongs.
 * derivId: ID of the unknown quantity with respect to which the
 *          equation should be differentiated.
 * jac:     Jacobian matrix (block) to set.
 * x:       Value of the unknown quantity.
 */
void Equation::SetJacobianBlockBC(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x
) {
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);
}

/**
 * Set the linear operator matrix elements in the given
 * matrix in order to represent this equation.
 *
 * mat: Matrix to set elements of.
 * rhs: Vector representing equation right-hand-side.
 */
void Equation::SetMatrixElements(Matrix *mat, real_t *rhs) {
    if (this->IsPredetermined()) {
        this->predetermined->SetMatrixElements(mat, rhs);
    } else {
        for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
            (*it)->SetMatrixElements(mat, rhs);

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
}

/**
 * Evaluate this equation and assign its value to the
 * given function vector.
 * 
 * vec: Function vector to assign evaluated equation to.
 * x:   Value of the unknown to evaluate the function for.
 */
void Equation::SetVectorElements(real_t *vec, const real_t *x) {
    if (this->IsPredetermined()) {
        this->predetermined->SetVectorElements(vec, x);
    } else {
        for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
            (*it)->SetVectorElements(vec, x);

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
}

