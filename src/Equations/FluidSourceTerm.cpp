/**
 * Implementation of a base class for source terms (kinetic or fluid) which 
 * are proportional to a fluid unknown quantity, of the form
 *     T = S(r,p1,p2,y) * x(r,t),
 * where x is the fluid unknown and S describes an arbitrary source
 * function which may be a (local) function of unknowns y.
 * Derived classes implement the evaluation of S and, if applicable,
 * the jacobian dS/dy. 
 * Is used for example by ParticleSourceTerm and AvalancheSourceRP.
 */

#include "DREAM/Equations/FluidSourceTerm.hpp"
#include <limits>

using namespace DREAM;

/**
 * Constructor
 */
FluidSourceTerm::FluidSourceTerm(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u
) : EquationTerm(kineticGrid), unknowns(u) {
    sourceVec = new real_t[kineticGrid->GetNCells()];
}

/**
 * Destructor
 */
FluidSourceTerm::~FluidSourceTerm(){
    delete [] sourceVec;
}

/**
 *  Reallocate when grid is rebuilt
 */
bool FluidSourceTerm::GridRebuilt(){
    delete [] sourceVec;
    sourceVec = new real_t[this->grid->GetNCells()];
    return true;    
}

/** 
 * Rebuild source vector
 */
void FluidSourceTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) {
    len_t offset=0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++)
                sourceVec[offset + j*n1[ir] + i] = GetSourceFunction(ir,i,j);
        offset += n1[ir]*n2[ir];
    }
}

/**
 * Normalizes the source vector so that it integrates 
 * over momentum to 'c' at each radius 
 */
void FluidSourceTerm::NormalizeSourceToConstant(const real_t c, real_t *normFactors){
    len_t offset=0;
    for(len_t ir = 0; ir<nr; ir++){
        real_t normFact = c/grid->IntegralMomentumAtRadius(ir,sourceVec+offset);
        len_t N = n1[ir]*n2[ir];
        for(len_t i=0; i<N; i++)
            sourceVec[offset+i] *= normFact;
        if(normFactors != nullptr)
            normFactors[ir] = normFact;
        offset += N; 
    }
}

/**
 * Set matrix elements.
 */
void FluidSourceTerm::SetMatrixElements(FVM::Matrix *mat, real_t* /*rhs*/){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                len_t ind = offset + n1[ir]*j + i;
                mat->SetElement(ind, ir, sourceVec[ind]);
            }
        offset += n1[ir]*n2[ir];
    }        
}

/**
 * Set vector elements.
 */
void FluidSourceTerm::SetVectorElements(real_t *vec, const real_t *x){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                len_t ind = offset + n1[ir]*j + i;
                vec[ind] += sourceVec[ind]*x[ir];
            }
        offset += n1[ir]*n2[ir];
    }
}

/**
 * Set jacobian matrix elements.
 */
bool FluidSourceTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x){
    bool contributes = (uqtyId == derivId);
    if(uqtyId == derivId)
        SetMatrixElements(jac, nullptr);

    // return if not derivId included in the list of dependent unknowns
    if(!HasJacobianContribution(derivId))
        return contributes;
    else
        contributes = true;

    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                real_t dS = GetSourceFunctionJacobian(ir,i,j,derivId);
                jac->SetElement(offset + n1[ir]*j + i, ir, dS*x[ir]);
            }
        offset += n1[ir]*n2[ir];
    }

    return contributes;
}
