/**
 * Implementation of a base class for equation terms that are a direct
 * weighted product of two unknown quantities:
 *      T = w * y * x,
 * where w is a grid-dependent weight and x and y are two unknown 
 * quantities defined on the same grid.
 * The input UnknownQuantityID wUqtyId refers to the quantity y,
 * and the created matrix is meant to "act on x".
 * Evaluation of weights must be implemented in derived classes. 
 */

#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/DiagonalQuadraticTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
DiagonalQuadraticTerm::DiagonalQuadraticTerm(Grid *g, const len_t uId, UnknownQuantityHandler *u)
        : DiagonalTerm(g) {
    this->wUqtyId = uId;
    this->unknowns = u;
}

/**
* if uqtyId = wUqty, this is done twice (once in SetMatrixElements 
* in SetJacobianBlock, and again here. Thus, quadratic terms of 
* the form x^2 correctly get a factor of 2. 
*/
void DiagonalQuadraticTerm::AddWeightsJacobian(
    const len_t /*uqtyId*/, const len_t derivId, Matrix *jac, const real_t* x
){
    if (derivId == wUqtyId){
        len_t N = this->grid->GetNCells();
        for (len_t i = 0; i < N; i++)
            jac->SetElement(i, i, x[i]*weights[i]);
    }
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void DiagonalQuadraticTerm::SetMatrixElements(Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();
    real_t *y = unknowns->GetUnknownData(wUqtyId);
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, y[i]*weights[i]);
}

/**
 * Set function vector for this term.
 */
void DiagonalQuadraticTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->grid->GetNCells();
    real_t *y = unknowns->GetUnknownData(wUqtyId);
    for (len_t i = 0; i < N; i++)
        vec[i] += weights[i] * y[i] * x[i];
}

