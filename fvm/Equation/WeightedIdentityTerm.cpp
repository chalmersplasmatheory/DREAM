/**
 * Implementation of a term which represents an identity term 
 * multiplied by an arbitrary grid-dependent weight function. 
 * Evaluation of weights must be implemented in derived classes. 
 */

#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/WeightedIdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
WeightedIdentityTerm::WeightedIdentityTerm(Grid *g)
        : EvaluableEquationTerm(g) { 
    //GridRebuilt();
//    weights = new real_t[grid->GetNCells()];
//    if(!TermDependsOnUnknowns())
//        SetWeights();
}

/**
 * Destructor
 */
WeightedIdentityTerm::~WeightedIdentityTerm() {
    this->DeallocateMemory();
    if(weights!=nullptr)
        delete [] weights;
    }

/**
 * Called if the grid is rebuilt; reallocates and rebuilds quantities.
 */
bool WeightedIdentityTerm::GridRebuilt(){
    this->AllocateMemory();
    if(weights!=nullptr)
        delete [] weights;
    weights = new real_t[grid->GetNCells()];
    if(!TermDependsOnUnknowns())
        SetWeights();
    return true;
}


/**
 * Evaluate this identity term.
 */
real_t WeightedIdentityTerm::Evaluate(real_t *vec, const real_t *x, const len_t eqnId, const len_t uqtyId) {
    if (eqnId == uqtyId) {
        throw EquationException(
            "The weighted identity term is not intended to be used in the equation "
            "for the unknown quantity to which it is applied."
        );
        //return 1;
    }

    /*len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        vec[i] -= weights[i]*x[i];*/
    this->SetVectorElements(vec, x);

    return 1;
}

/**
 * Set a block for this term in the given jacobian matrix.
 */
void WeightedIdentityTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* /*x*/
) {
    if (derivId == uqtyId) {
        this->SetMatrixElements(jac, nullptr);
    }
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void WeightedIdentityTerm::SetMatrixElements(Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, weights[i]);
}

/**
 * Set function vector for this term.
 */
void WeightedIdentityTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        vec[i] += weights[i] * x[i];
}

