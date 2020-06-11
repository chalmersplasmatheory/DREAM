/**
 * Implementation of a base class for equation terms that are a direct
 * weighted product of two unknown quantities:
 *      T = w * y * x,
 * where w is a grid-dependent weight and x and y are two unknown 
 * quantities defined on the same grid.
 * The input UnknownQuantityID wUqtyId refers to the quantity y,
 * and the created matrix is meant to "act on x".
 * Evaluation of weights must be implemented in derived classes. 
 * The secondary quantity wUqtyId can be ions, in which case the
 * number of charge states are kept in wUqtyNMultiples. In that case,
 * weights are assumed to be 
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
    this->unknowns = u;
    this->wUqtyId = uId;
    this->wUqtyNMultiples = unknowns->GetUnknown(wUqtyId)->NumberOfMultiples();
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
            for(len_t n=0; n<wUqtyNMultiples; n++)            
                jac->SetElement(i, n*N+i, x[i]*weights[n*N+i]);
    }
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void DiagonalQuadraticTerm::SetMatrixElements(Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();
    real_t *y = unknowns->GetUnknownData(wUqtyId);

    for (len_t i = 0; i < N; i++)
        for(len_t n=0; n<wUqtyNMultiples; n++)
            mat->SetElement(i, i, y[N*n+i]*weights[N*n+i]);
}

/**
 * Set function vector for this term.
 */
void DiagonalQuadraticTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->grid->GetNCells();
    real_t *y = unknowns->GetUnknownData(wUqtyId);
    for (len_t i = 0; i < N; i++)
        for(len_t n=0; n<wUqtyNMultiples; n++)
            vec[i] += weights[n*N+i] * y[n*N+i] * x[i];
}

