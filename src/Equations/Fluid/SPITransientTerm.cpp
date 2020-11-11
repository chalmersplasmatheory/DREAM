/**
 * Implementation of an Euler backward transient term multiplied
 * by an arbitrary grid-dependent weight function. 
 * Evaluation of weights must be implemented in derived classes. 
 */

#include <iostream>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "DREAM/Equations/Fluid/SPITransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
SPITransientTerm::SPITransientTerm(FVM::Grid *grid, const len_t unknownId, len_t nShard, real_t scaleFactor)
        : FVM::EquationTerm(grid), unknownId(unknownId), nShard(nShard), scaleFactor(scaleFactor){}

/**
 * Rebuild the transient term.
 *
 * dt: Length of next time step to take.
 */
void SPITransientTerm::Rebuild(const real_t, const real_t dt, FVM::UnknownQuantityHandler *uqty) {
    this->dt = dt;
    this->xn = uqty->GetUnknownDataPrevious(this->unknownId);
}

void SPITransientTerm::SetJacobianBlock(const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(derivId==this->unknownId){
        for (len_t i = 0; i < nShard; i++){
            //jac->SetElement(i,i, 5.0/3.0*scaleFactor*pow(xn[i],2.0/3.0)/this->dt);
            jac->SetElement(i,i, 5.0/9.0*scaleFactor*pow(xn[i],-4.0/9.0)/this->dt);
        }
    }
}

void SPITransientTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    this->SetVectorElements(rhs, nullptr);
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * xnp1:   Previous solution.
 */
void SPITransientTerm::SetVectorElements(real_t *vec, const real_t *xnp1) {

    for (len_t i = 0; i < nShard; i++)
        //vec[i] += scaleFactor*(pow(xnp1[i],5.0/3.0) - pow(xn[i],5.0/3.0)) / this->dt;
        vec[i] += scaleFactor*(pow(xnp1[i],5.0/9.0) - pow(xn[i],5.0/9.0)) / this->dt;
}

