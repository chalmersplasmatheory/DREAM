/**
 * Implementation of a base class for equation terms with diagonal
 * matrices. 
 * Evaluation of weights must be implemented in derived classes. 
 */

#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/DiagonalTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
DiagonalTerm::DiagonalTerm(Grid *g) : EquationTerm(g){
//    InitializeWeights();
}

/**
 * Destructor
 */
DiagonalTerm::~DiagonalTerm() {
    this->DeallocateWeights();
}


/**
 * Allocate and set weights.
 */
void DiagonalTerm::InitializeWeights(){
    AllocateWeights(); 
    SetWeights();
}

/**
 * If term depends on unknowns, set weights.
 * Otherwise this is done only in GridRebuilt.
 */
void DiagonalTerm::Rebuild(const real_t, const real_t, UnknownQuantityHandler*){ 
    if(!hasBeenInitialized){
        InitializeWeights();
        hasBeenInitialized = true;
    }
    if(TermDependsOnUnknowns()) 
        SetWeights();
}


/**
 * Called if the grid is rebuilt; reallocates and rebuilds quantities.
 */
bool DiagonalTerm::GridRebuilt(){
    this->AllocateMemory();
    this->AllocateWeights();
    if(!TermDependsOnUnknowns())
        SetWeights();
    return true;
}


/**
 * Evaluate this identity term.
 */
/*
real_t DiagonalTerm::Evaluate(real_t *vec, const real_t *x, const len_t eqnId, const len_t uqtyId) {
    if (eqnId == uqtyId) {
        throw EquationException(
            "The diagonal term is not intended to be used in the equation "
            "for the unknown quantity to which it is applied."
        );
        //return 1;
    }

    this->SetVectorElements(vec, x);

    return 1;
}
*/
/**
 * Set a block for this term in the given jacobian matrix.
 */
void DiagonalTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* x
) {
    if (derivId == uqtyId) {
        this->SetMatrixElements(jac, nullptr);
    }
    AddWeightsJacobian(uqtyId, derivId, jac, x);
}

/**
 * Allocate weights
 */
void DiagonalTerm::AllocateWeights(){
    DeallocateWeights(); 
    weights = new real_t[grid->GetNCells()];
}

/**
 * Deallocate weights
 */
void DiagonalTerm::DeallocateWeights(){
    if(weights!=nullptr) 
    delete[] weights;
}


