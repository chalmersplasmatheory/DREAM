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
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
DiagonalComplexTerm::DiagonalComplexTerm(Grid *g, UnknownQuantityHandler *u)
        : DiagonalTerm(g) {
    this->unknowns = u;
    AllocateDiffWeights();
}

DiagonalComplexTerm::~DiagonalComplexTerm(){
    DeallocateDiffWeights();
}

/**
* if uqtyId = wUqty, this is done twice (once in SetMatrixElements 
* in SetJacobianBlock, and again here. Thus, quadratic terms of 
* the form x^2 correctly get a factor of 2. 
*/
void DiagonalComplexTerm::AddWeightsJacobian(
    const len_t /*uqtyId*/, const len_t derivId, Matrix *jac, const real_t* x
){
    /**
    * Check if derivId is one of the id's that contributes 
    * to this advection coefficient 
    */
    bool hasDerivIdContribution = false;
    len_t nMultiples;
    for(len_t i_deriv = 0; i_deriv < derivIds.size(); i_deriv++){
        if (derivId == derivIds[i_deriv]){
            nMultiples = derivNMultiples[i_deriv];
            hasDerivIdContribution = true;
        }
    }
    if(!hasDerivIdContribution)
        return;
    
    SetDiffWeights(derivId, nMultiples);

    len_t NCells = grid->GetNCells();
    for(len_t n=0; n<nMultiples; n++)
        for(len_t i=0; i<NCells; i++)
            jac->SetElement(i, n*NCells+i, diffWeights[n*NCells + i] * x[i] ); 

    ResetDiffWeights();

}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void DiagonalComplexTerm::SetMatrixElements(Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, weights[i]);
}

/**
 * Set function vector for this term.
 */
void DiagonalComplexTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        vec[i] += weights[i] * x[i];
}


/**
 * Set all diffweights to 0.
 */
void DiagonalComplexTerm::ResetDiffWeights(){
    len_t nMultiples = MaxNMultiple();
    len_t NCells = grid->GetNCells();

    for(len_t i = 0; i<NCells*nMultiples; i++){
        diffWeights[i] = 0;
    }
}

/**
 * Allocate differentiation coefficients.
 */
void DiagonalComplexTerm::AllocateDiffWeights() {
    DeallocateDiffWeights();
    len_t nMultiples = MaxNMultiple();
    len_t NCells = grid->GetNCells();
    diffWeights = new real_t[NCells*nMultiples];
    ResetDiffWeights();
}

/**
 * Deallocate differentiation coefficients.
 */
void DiagonalComplexTerm::DeallocateDiffWeights() {
    delete [] diffWeights;
}

