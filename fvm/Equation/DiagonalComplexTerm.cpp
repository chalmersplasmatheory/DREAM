/**
 * Implementation of a base class for diagonal equation terms 
 * where the coefficients are a complex function of unknowns:
 *      T = w(U) * x,
 * where w are the weights and U represents all unknown quantities. 
 * Evaluation of weights must be implemented in derived classes. 
 * It is essentially a DiagonalLinearTerm, but where extra support
 * is provided for setting the Jacobian via SetDiffWeights.
 */

#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
DiagonalComplexTerm::DiagonalComplexTerm(Grid *g, UnknownQuantityHandler *u)
        : DiagonalTerm(g) {
    this->unknowns = u;
}

DiagonalComplexTerm::~DiagonalComplexTerm(){
    DeallocateDiffWeights();
}

/**
* SetDiffWeights sets the weight jacobian dw/dUderiv, assumed to be local
* (so that w at phase-space point z only depends on U(z) and not, for example,
* on integrals of U). 
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
    for(len_t i_deriv = 0; i_deriv < derivIds.size(); i_deriv++)
        if (derivId == derivIds[i_deriv]){
            nMultiples = derivNMultiples[i_deriv];
            hasDerivIdContribution = true;
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
    if(diffWeights != nullptr)
        delete [] diffWeights;
}

