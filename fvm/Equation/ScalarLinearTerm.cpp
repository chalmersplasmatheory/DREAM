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
    : EquationTerm(scalarGrid), uqtyId(uqtyId){
    this->targetGrid = targetGrid;
    this->unknowns = u;
    AllocateWeights();
}

/**
 * Destructor.
 */
ScalarLinearTerm::~ScalarLinearTerm() {
	DeallocateWeights();
}

/**
 * Transform the given input vector 'vec' so as to solve for
 * the value of the unknown quantity operated on by this
 * ScalarLinearTerm. We imagine that this term is part of
 * an equation of the form
 *
 *   w*x + F(y) = 0
 *
 * where F denotes an arbitrary, evaluable function of any
 * other unknown quantities (i.e. y may NOT be x). The input
 * vector is assumed to contain the value 'F', and we
 * subsequently transform the input vector so that the output
 * represents
 *
 *   vec[i] = -F_i(y) / w_i
 *
 * i.e. the corresponding value of the unknown quantity 'x'.
 *
 * vec: Vector, containing the values corresponding to 'F'
 *      to transform.
 */
/*
void ScalarLinearTerm::EvaluableTransform(real_t *vec) {
    const len_t N = this->grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] = -vec[i] / weights[i];
}
*/

/**
 * Rebuild the weights for this scalar linear term.
 */
void ScalarLinearTerm::Rebuild(const real_t, const real_t, UnknownQuantityHandler*){
    SetWeights();
}

/**
 * Method called whenever the associated 'scalarGrid' is
 * rebuilt.
 */
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
bool ScalarLinearTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* /*x*/
) {
    bool contributes = false;
    if (derivId == uqtyId) {
        this->SetMatrixElements(jac, nullptr);
        contributes = true;
    }

    return contributes;
}

/**
 * Allocate memory for weights
 */
void ScalarLinearTerm::AllocateWeights(){
    DeallocateWeights();
    nWeights = targetGrid->GetNCells() * unknowns->GetUnknown(uqtyId)->NumberOfMultiples();
    weights = new real_t[nWeights];
    for(len_t i=0; i<nWeights; i++)
        weights[i] = 0;
}

/**
 * Deallocate memory for weights
 */
void ScalarLinearTerm::DeallocateWeights(){
    if(weights != nullptr)
        delete [] weights;
}
