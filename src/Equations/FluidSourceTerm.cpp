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
) : EquationTerm(kineticGrid), unknowns(u) {}

/**
 * Set matrix elements.
 */
void FluidSourceTerm::SetMatrixElements(FVM::Matrix *mat, real_t* /*rhs*/){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                real_t S = GetSourceFunction(ir,i,j);
                mat->SetElement(offset + n1[ir]*j + i, ir, S);
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
                real_t S = GetSourceFunction(ir,i,j);
                vec[offset + n1[ir]*j + i] += S*x[ir];
            }
        offset += n1[ir]*n2[ir];
    }
}

/**
 * Set jacobian matrix elements.
 */
void FluidSourceTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x){
    if(uqtyId == derivId)
        SetMatrixElements(jac, nullptr);

    // check whether derivId is included in the list of dependent unknowns
    bool hasDerivIdContribution = false;
    for(len_t i_deriv = 0; i_deriv < derivIds.size(); i_deriv++)
        if (derivId == derivIds[i_deriv])
            hasDerivIdContribution = true;
    // if not: return
    if(!hasDerivIdContribution)
        return;


    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                real_t dS = GetSourceFunctionJacobian(ir,i,j,derivId);
                jac->SetElement(offset + n1[ir]*j + i, ir, dS*x[ir]);
            }
    offset += n1[ir]*n2[ir];
    }   
}
