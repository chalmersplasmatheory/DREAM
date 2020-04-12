/**
 * Implementation of a general moment function
 * (such as the moment of a distribution function).
 */

#include "FVM/Equation/MomentQuantity.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
MomentQuantity::MomentQuantity(Grid *momentGrid, Grid *fGrid, len_t momentId, len_t fId) 
    : EquationTerm(momentGrid), fGrid(fGrid), momentId(momentId), fId(fId) {
    
    this->integrand = new real_t[fGrid->GetNCells()];

    // Figure out the maximum number of non-zeros needed per
    // matrix row...
    const len_t nr = fGrid->GetNr();
    this->nnz_per_row = 0;
    for (len_t i = 0; i < nr; i++) {
        len_t nc = fGrid->GetMomentumGrid(i)->GetNCells();
        if (this->nnz_per_row < nc)
            this->nnz_per_row = nc;
    }
}

/**
 * Destructor.
 */
MomentQuantity::~MomentQuantity() {
    delete [] this->integrand;
}


/**
 * Set the jacobian elements for this term. This assumes that
 * the integrand does not depend on any unknown quantity (other than
 * linearly on the unknown to which this moment operator is applied
 * (usually the distribution function f)).
 *
 * derivId: Unknown ID of derivative with respect to which differentiation
 *          should be done.
 * unknId:  ID of the unknown to differentiate.
 * jac:     Jacobian matrix to set elements of.
 */
void MomentQuantity::SetJacobianBlock(
    const len_t derivId, const len_t unknId, Matrix *jac
) {
    if (derivId == fId && unknId == momentId) {
        #define X(IR,I,J,V) jac->SetElement(offset+((J)*np1+(I)), offset+((J)*np1+(I)), (V))
        #   include "MomentQuantity.setel.cpp"
        #undef X
    }
}

/**
 * Set the elements of the linear operator matrix corresponding to
 * this operator.
 *
 * mat: Linear operator matrix to set elements of.
 * rhs: Equation right-hand-side.
 */
void MomentQuantity::SetMatrixElements(Matrix *mat, real_t*) {
    #define X(IR,I,J,V) mat->SetElement((IR), offset + ((J)*np1 + (I)), (V))
    #   include "MomentQuantity.setel.cpp"
    #undef X
}

/**
 * Set the elements of the function vector 'F' in the non-linear
 * solver.
 *
 * vec: Vector to set elements of.
 * f:   Current value of the unknown quantity for which to evaluate
 *      this operator.
 */
void MomentQuantity::SetVectorElements(real_t *vec, const real_t *f) {
    #define X(IR,I,J,V) vec[(IR)] += f[offset+((J)*np1+(I))] * (V)
    #   include "MomentQuantity.setel.cpp"
    #undef X
}

