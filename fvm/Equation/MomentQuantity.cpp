/**
 * Implementation of a general moment function
 * (such as the moment of a distribution function).
 */

#include "FVM/Equation/MomentQuantity.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
MomentQuantity::MomentQuantity(Grid *momentGrid, Grid *fGrid, len_t momentId, len_t fId, UnknownQuantityHandler *u, real_t pThreshold, pThresholdMode pMode) 
    : EquationTerm(momentGrid), fGrid(fGrid), momentId(momentId), 
      fId(fId), unknowns(u), pThreshold(pThreshold), pMode(pMode) {
    this->hasThreshold = (pThreshold!=0);
    this->GridRebuilt();
}

/**
 * Destructor.
 */
MomentQuantity::~MomentQuantity() {
    delete [] this->integrand;
}

/**
 * Method that is called whenever the grid is rebuilt.
 * Here, we use it to rebuild the 'integrand' variable.
 */
bool MomentQuantity::GridRebuilt() {
    bool rebuilt = this->EquationTerm::GridRebuilt();
    
    const len_t N = this->fGrid->GetNCells();

    if (this->nIntegrand != N) {
        this->nIntegrand = N;
        this->integrand = new real_t[N];

        // Figure out the maximum number of non-zeros needed per
        // matrix row...
        const len_t nr = fGrid->GetNr();
        this->nnz_per_row = 0;
        for (len_t i = 0; i < nr; i++) {
            len_t nc = fGrid->GetMomentumGrid(i)->GetNCells();
            if (this->nnz_per_row < nc)
                this->nnz_per_row = nc;
        }

        return true;
    } else return rebuilt;
}

/**
 * Tests whether a given momentum satisfies the provided thresholds.
 * If not, the integrand will be set to zero in those points.
 * If a non-zero threshold is set and MIN limit, at least one momentum
 * grid point will always be assumed to contribute.
 * Modes:
 *  MC: assumes pThreshold was given in units of m*c
 *  THERMAL: assumes pThreshold was given in thermal electron momenta
 *  MIN: assumes pThreshold is a lower limit
 *  MAX: assumes pThreshold is an upper limit
 */
bool MomentQuantity::SatisfiesThreshold(len_t ir, len_t i1, len_t i2){
    if(!this->hasThreshold)
        return true;

    const real_t p = fGrid->GetMomentumGrid(ir)->GetP(i1,i2);
    switch(pMode){
        case P_THRESHOLD_MODE_MIN_MC: // XXX: assumes p-xi grid
            return (p>pThreshold) || ((i1==0) && (pThreshold!=0));
        case P_THRESHOLD_MODE_MAX_MC:
            return p<pThreshold;
        case P_THRESHOLD_MODE_MIN_THERMAL:{
            const real_t Tcold = unknowns->GetUnknownData(OptionConstants::UQTY_T_COLD)[ir];
            return (p > pThreshold * sqrt(2*Tcold/Constants::mc2inEV))  || ((i1==0) && (pThreshold!=0));
        }
        case P_THRESHOLD_MODE_MAX_THERMAL:{
            const real_t Tcold = unknowns->GetUnknownData(OptionConstants::UQTY_T_COLD)[ir];
            return p < pThreshold * sqrt(2*Tcold/Constants::mc2inEV);
        }
        default:
            throw FVM::FVMException("MomentQuantity: Unrecognized p threshold mode.");
            return 0;
    }

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
 * x:       Value of the unknown quantity.
 */
void MomentQuantity::SetJacobianBlock(
    const len_t unknId, const len_t derivId, Matrix *jac, const real_t* /*x*/
) {
    if (derivId == fId && unknId == fId)
        this->SetMatrixElements(jac,nullptr);
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

