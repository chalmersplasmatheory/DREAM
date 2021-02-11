/**
 * Implementation of a general moment function
 * (such as the moment of a distribution function).
 */

#include "FVM/Equation/MomentQuantity.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

#include <iostream>

using namespace DREAM::FVM;


/**
 * Constructor.
 */
MomentQuantity::MomentQuantity(Grid *momentGrid, Grid *fGrid, len_t momentId, len_t fId, UnknownQuantityHandler *u, real_t pThreshold, pThresholdMode pMode) 
    : EquationTerm(momentGrid), fGrid(fGrid), momentId(momentId), 
      fId(fId), unknowns(u), pThreshold(pThreshold), pMode(pMode) {
    this->hasThreshold = (pThreshold!=0);
    id_Tcold = u->GetUnknownID(OptionConstants::UQTY_T_COLD);
    if ((pMode == P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH) || (pMode == P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH)){
        AddUnknownForJacobian(u, id_Tcold);
        smoothEnvelopeStepWidth = 2;
    }
    this->GridRebuilt();    
}

/**
 * Destructor.
 */
MomentQuantity::~MomentQuantity() {
    delete [] this->integrand;
    delete [] this->diffIntegrand;
}

/**
 * Rebuild allocates memory for the integrand and diffIntegrand
 * arrays whenever the grid has changed size.
 */
bool MomentQuantity::GridRebuilt() {
    bool rebuilt = this->EquationTerm::GridRebuilt();
    const len_t N = this->fGrid->GetNCells();
    if (this->nIntegrand != N) {
        this->nIntegrand = N;
        if(integrand != nullptr)
            delete [] integrand;
        this->integrand = new real_t[N];

        if(GetMaxNumberOfMultiplesJacobian())
            AllocateDiffIntegrand();

        return true;
    } else
        return rebuilt;   
}

void MomentQuantity::AllocateDiffIntegrand(){
    len_t nMultiples = GetMaxNumberOfMultiplesJacobian();
    if(nMultiples){
        if(diffIntegrand != nullptr)
            delete [] diffIntegrand;    
        const len_t N = this->fGrid->GetNCells();
        this->diffIntegrand = new real_t[nMultiples*N];
    }
}



/**
 * returns the momentum grid spacing dp at momentum p0
 */
real_t FindThresholdStep(real_t p0, MomentumGrid *mg){
    const real_t *p = mg->GetP1();
    for(len_t i=0; i<mg->GetNp1(); i++)
        if(p[i]>p0)
            return mg->GetDp1(i);
    return p0 - mg->GetP1(mg->GetNp1()-1); // what else to do if p0 is outside the grid? 
}

/**
 * An envelope added to the integrand, which can be used to set thresholds
 * and integration limits, or more general modifications.
 * Modes:
 *  MC: assumes pThreshold was given in units of m*c
 *  THERMAL: assumes pThreshold was given in thermal electron momenta
 *  MIN: assumes pThreshold is a lower limit
 *  MAX: assumes pThreshold is an upper limit
 *  SMOOTH: changes the limit to a smooth tanh step, with a width of 
 *          smoothEnvelopeStepWidth grid points in each direction around pThreshold
 * 
 * XXX: Assumes p-xi grid
 */
real_t MomentQuantity::ThresholdEnvelope(len_t ir, len_t i){
    if(!this->hasThreshold)
        return 1;
    MomentumGrid *mg = fGrid->GetMomentumGrid(ir);

    const real_t 
        p   = mg->GetP1(i),
        p_u = mg->GetP1_f(i+1),
        p_l = mg->GetP1_f(i);
    
    // fracCellInRegion is 0 outside the region, 1 inside the 
    // region and the fraction of overlap dpOverlap/dp when
    // the cell straddles the region boundary. 
    switch(pMode){
        case P_THRESHOLD_MODE_MIN_MC:{
            real_t fracCellInRegion = 0;
            if(p_l>=pThreshold)
                fracCellInRegion=1;
            else if(p_u>=pThreshold)
                fracCellInRegion = (p_u-pThreshold)/(p_u-p_l);
            return fracCellInRegion;
        }
        case P_THRESHOLD_MODE_MIN_THERMAL:{
            const real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir];
            real_t p0 = pThreshold * sqrt(2*Tcold/Constants::mc2inEV);
            real_t fracCellInRegion = 0;
            if(p_l>=p0)
                fracCellInRegion=1;
            else if(p_u>=p0)
                fracCellInRegion = (p_u-p0)/(p_u-p_l);
            return fracCellInRegion;
        }
        case P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH:{
            const real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir];
            real_t p0 = pThreshold * sqrt(2*Tcold/Constants::mc2inEV);
            real_t dp = FindThresholdStep(p0, fGrid->GetMomentumGrid(ir));
            real_t x = (p-p0)/(smoothEnvelopeStepWidth*dp);
            real_t thx = std::tanh(x);
            return .5*( 1 + thx );
        }
        case P_THRESHOLD_MODE_MAX_MC:{
            real_t fracCellInRegion = 0;
            if(p_u<=pThreshold)
                fracCellInRegion=1;
            else if(p_l<pThreshold)
                fracCellInRegion = (pThreshold-p_l)/(p_u-p_l);
            return fracCellInRegion;
        }
        case P_THRESHOLD_MODE_MAX_THERMAL:{
            real_t fracCellInRegion = 0;
            const real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir];
            real_t p0 = pThreshold * sqrt(2*Tcold/Constants::mc2inEV); 
            if(p_u<=p0)
                fracCellInRegion=1;
            else if(p_l<p0)
                fracCellInRegion = (p0-p_l)/(p_u-p_l);
            return fracCellInRegion;
        }
        case P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH:{
            const real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir];
            real_t p0 = pThreshold * sqrt(2*Tcold/Constants::mc2inEV);
            real_t dp = FindThresholdStep(p0, fGrid->GetMomentumGrid(ir));
            real_t x = (p-p0)/(smoothEnvelopeStepWidth*dp);
            real_t thx = std::tanh(-x);
            return .5*( 1 + thx );
        }
        default:
            throw FVM::FVMException("MomentQuantity: Unrecognized p threshold mode.");
            return 0;
    }
}

/**
 * Returns the jacobian with respect to Tcold[ir] of the smooth threshold functions
 * XXX: Assumes p-xi grids
 */
real_t MomentQuantity::DiffThresholdEnvelope(len_t ir, len_t i){
    if(!this->hasThreshold)
        return 0;
    MomentumGrid *mg = fGrid->GetMomentumGrid(ir);
    const real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir];        
    const real_t p = mg->GetP1(i);
    switch(pMode){
        case P_THRESHOLD_MODE_MIN_THERMAL:{
            const real_t pTe = sqrt(2*Tcold/Constants::mc2inEV);
            real_t p0 = pThreshold * pTe;
            real_t dp0 = pThreshold/(Constants::mc2inEV * pTe); 
            const real_t   
                p_u = mg->GetP1_f(i+1),
                p_l = mg->GetP1_f(i);
            real_t fracCellInRegion = 0;
            if(p_l<p0 && p_u>=p0)
                fracCellInRegion = -dp0/(p_u-p_l);
            return fracCellInRegion;
        }
        case P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH:{
            // XXX: assumes p-xi grid
            real_t p0 = pThreshold * sqrt(2*Tcold/Constants::mc2inEV);
            real_t dp = FindThresholdStep(p0, mg);
            real_t x = (p-p0)/(smoothEnvelopeStepWidth*dp);
            real_t chx = std::cosh(x);
            if(chx>sqrt(std::numeric_limits<real_t>::max()))
                return 0;
            else
                return -p0/(4*Tcold*smoothEnvelopeStepWidth*dp*chx*chx);
        }
        case P_THRESHOLD_MODE_MAX_THERMAL:{
            const real_t pTe = sqrt(2*Tcold/Constants::mc2inEV);
            real_t p0 = pThreshold * pTe;
            real_t dp0 = pThreshold/(Constants::mc2inEV * pTe); 
            const real_t   
                p_u = mg->GetP1_f(i+1),
                p_l = mg->GetP1_f(i);
            real_t fracCellInRegion = 0;
            if(p_l<p0 && p_u>=p0)
                fracCellInRegion = dp0/(p_u-p_l);
            return fracCellInRegion;
        }
        case P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH:{
            // XXX: assumes p-xi grid
            const real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir];
            real_t p0 = pThreshold * sqrt(2*Tcold/Constants::mc2inEV);
            real_t dp = FindThresholdStep(p0, mg);
            real_t x = (p-p0)/(smoothEnvelopeStepWidth*dp);
            real_t chx = std::cosh(x);
            if(chx>sqrt(std::numeric_limits<real_t>::max()))
                return 0;
            else
                return p0/(4*Tcold*smoothEnvelopeStepWidth*dp*chx*chx);
        }
        default:
            return 0;
    }

}

/**
 * Adds the jacobian of the envelope to the integrand jacobian.
 */
void MomentQuantity::AddDiffEnvelope(){
    len_t offset=0;
    for(len_t ir=0; ir<fGrid->GetNr(); ir++){
        MomentumGrid *mg = fGrid->GetMomentumGrid(ir);
        len_t np1 = mg->GetNp1();
        len_t np2 = mg->GetNp2();
        for(len_t i=0;i<np1;i++)
            for(len_t j=0; j<np2;j++){
                len_t idx = np2*j + i;
                diffIntegrand[offset + idx] += 
                    (DiffThresholdEnvelope(ir,i)/ThresholdEnvelope(ir,i)) 
                    * integrand[offset + idx];
            }
        offset += np1*np2;
    }
}

/**
 * Set the jacobian elements for this term.
 *
 * derivId: Unknown ID of derivative with respect to which differentiation
 *          should be done.
 * unknId:  ID of the unknown to differentiate.
 * jac:     Jacobian matrix to set elements of.
 * x:       Value of the unknown quantity.
 */
void MomentQuantity::SetJacobianBlock(
    const len_t unknId, const len_t derivId, Matrix *jac, const real_t* x
) {
    if (derivId == unknId)
        this->SetMatrixElements(jac,nullptr);

    // Add potential contributions from fluid quantities 
    bool hasDerivIdContribution = false;
    len_t nMultiples;
    for(len_t i_deriv = 0; i_deriv < derivIds.size(); i_deriv++)
        if (derivId == derivIds[i_deriv]){
            nMultiples = derivNMultiples[i_deriv];
            hasDerivIdContribution = true;
        }

    if(!hasDerivIdContribution)
        return;

    SetDiffIntegrand(derivId);
    if((derivId==id_Tcold) && ((pMode == P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH) || (pMode == P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH)))
        AddDiffEnvelope();
    
    len_t offset_n = 0;
//    #define X(IR,I,J,V) jac->SetElement((IR), (IR)+n*nr, (V)*x[offset+((J)*np1+(I))])
    #define X(ID,V) \
        VAL += (V)*x[offset+(ID)]*diffIntegrand[offset_n + (ID)];
    #define ApplyX(IR) \
        IND = (IR) + n*nr; \
        jac->SetRow((IR), 1, &IND, &VAL); \
        VAL = 0;
    PetscInt IND;
    PetscScalar VAL = 0;
    for(len_t n=0; n<nMultiples;n++){
        #include "MomentQuantity.setel.cpp"
        offset_n += offset; 
    }
    #undef ApplyX
    #undef X
}

/**
 * Set the elements of the linear operator matrix corresponding to
 * this operator.
 *
 * mat: Linear operator matrix to set elements of.
 * rhs: Equation right-hand-side.
 */
void MomentQuantity::SetMatrixElements(Matrix *mat, real_t*) {
    #define X(ID,V) \
        IND_ARR[N_IND] = offset + (ID); \
        VAL_ARR[N_IND]= (V) * integrand[offset + (ID)]; \
        N_IND++; 
    #define ApplyX(IR) \
        mat->SetRow((IR),N_IND,IND_ARR, VAL_ARR); \
        N_IND = 0;
    len_t N_IND = 0;
    const len_t SIZE = GetNumberOfNonZerosPerRow();
    PetscInt *IND_ARR = new PetscInt[SIZE];
    PetscScalar *VAL_ARR = new PetscScalar[SIZE];
    #   include "MomentQuantity.setel.cpp"
    delete [] IND_ARR;
    delete [] VAL_ARR;
    #undef ApplyX
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
    #define X(ID,V) vec[ir] += f[offset+(ID)] * integrand[offset+(ID)] * (V);
    #define ApplyX(IR) {}; // do nothing
    #   include "MomentQuantity.setel.cpp"
    #undef ApplyX
    #undef X
}

