/**
 * Implementation of a term which represents F(r,p1,p2) times the quantity 
 * it is applied to, where F is described by the lambda function 
 * weightFunc = [](len_t ir, len_t i, len_t j){return F(r_ir,p1_i,p2_j)}
 */

#include "FVM/Equation/WeightedIdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
WeightedIdentityTerm::WeightedIdentityTerm(Grid *g, std::function<real_t(len_t,len_t,len_t)> *weightFunc)
        : EvaluableEquationTerm(g), weightFunc(weightFunc) { 
    weights = new real_t[grid->GetNCells()];
}

/**
 * Destructor
 */
WeightedIdentityTerm::~WeightedIdentityTerm() {
    delete [] weights;
}

void WeightedIdentityTerm::Rebuild(const real_t, const real_t, UnknownQuantityHandler*){
    len_t offset = 0;
    std::function<real_t(len_t,len_t,len_t)> func = *weightFunc;
    for (len_t ir = 0; ir < nr; ir++)
        for(len_t i = 0; i < n1[ir]; i++)
            for(len_t j = 0; j < n2[ir]; j++){
                weights[offset + n1[ir]*j + i] = func(ir,i,j);
                offset += n1[ir]*n2[ir];
            }

}



/**
 * Evaluate this identity term.
 */
void WeightedIdentityTerm::Evaluate(real_t *vec, const real_t *x, const len_t eqnId, const len_t uqtyId) {
    if (eqnId == uqtyId)
        return;

    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        vec[i] -= weights[i]*x[i];
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
        vec[i] = weights[i] * x[i];
}

