/**
 * Implementation of a base class for scalar equation terms that are linear,
 * T = W*x = sum_i w_i x_i
 * where x is an unknown quantity and W is a constant operator 
 * (which can be represented by a row matrix of same length as x) 
 */

#include "FVM/Equation/ScalarLinearTerm.hpp"


using namespace DREAM::FVM;

/**
 * Constructor. targetGrid is the grid of the unknown quantity
 * of which the term "acts on", with uqtyId its ID. 
 */
ScalarLinearTerm::ScalarLinearTerm(Grid *scalarGrid,Grid *targetGrid,
    UnknownQuantityHandler *u, const len_t uqtyId)
    : EvaluableEquationTerm(scalarGrid), uqtyId(uqtyId){
    this->targetGrid = targetGrid;
    this->unknowns = u;
    AllocateWeights();
}

void ScalarLinearTerm::Rebuild(const real_t, const real_t, UnknownQuantityHandler*){
    SetWeights();
}

bool ScalarLinearTerm::GridRebuilt(){
    AllocateMemory();
    AllocateWeights();

    return true;
}


/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void ScalarLinearTerm::SetMatrixElements(Matrix *mat, real_t*) {
    for (len_t i = 0; i < nWeights; i++)
        mat->SetElement(0, i, weights[i]);
}

/**
 * Set function vector for this term.
 */
void ScalarLinearTerm::SetVectorElements(real_t *vec, const real_t *x) {
    for (len_t i = 0; i < nWeights; i++)
        vec[0] += weights[i] * x[i];
}

/**
 * Set a block for this term in the given jacobian matrix.
 */
void ScalarLinearTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* /*x*/
) {
    if (derivId == uqtyId) {
        this->SetMatrixElements(jac, nullptr);
    }
}

real_t* ScalarLinearTerm::Evaluate(real_t *vec, const real_t *x, const len_t eqnId, const len_t uqtyId) {
    // The diagonal identity term is implied when evaluating equation terms,
    // and so we shouldn't add it here...
    if ((eqnId == uqtyId) && (nWeights==1))
        return weights;

    /*len_t N = this->grid->GetNCells();
    const real_t sf = this->scaleFactor;
    for (len_t i = 0; i < N; i++)
        vec[i] -= sf*x[i];*/
    this->SetVectorElements(vec, x);

    return nullptr;
}




void ScalarLinearTerm::AllocateWeights(){
    DeallocateWeights();
    nWeights = targetGrid->GetNCells() * unknowns->GetUnknown(uqtyId)->NumberOfMultiples();
    weights = new real_t[nWeights];
    for(len_t i=0; i<nWeights; i++)
        weights[i] = 0;
}

void ScalarLinearTerm::DeallocateWeights(){
    if(weights != nullptr)
        delete [] weights;
}
