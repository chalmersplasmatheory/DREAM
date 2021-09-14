/**
 * Implementation of the transient term for the SPI shard radii.
 * We make special implementation for this transient term since the shard radii variable has many multiples
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

bool SPITransientTerm::SetJacobianBlock(const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(derivId==this->unknownId){
        for (len_t i = 0; i < nShard; i++){
            jac->SetElement(i,i, scaleFactor/this->dt);
        }
        return true;
    }
    return false;
}

void SPITransientTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    this->SetVectorElements(rhs, nullptr);
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * xnp1:   Current solution.
 */
void SPITransientTerm::SetVectorElements(real_t *vec, const real_t *xnp1) {

    for (len_t i = 0; i < nShard; i++)
        vec[i] += scaleFactor*(xnp1[i] - xn[i]) / this->dt;
}

