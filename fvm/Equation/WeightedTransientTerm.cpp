/**
 * Implementation of a Euler backward transient term with weight function F.
 *
 *     df       f_{n+1} - f_n
 * F * -- ~ F * -------------
 *     dt            dt
 * where F = F(r,p1,p2).
 * F is described by the lambda function weightFunc = [](len_t ir, len_t i, len_t j){return F(r_ir,p1_i,p2_j)}
 */

#include <iostream>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Equation/WeightedTransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM::FVM;


/**
 * Constructor.
 */
WeightedTransientTerm::WeightedTransientTerm(Grid *grid, const len_t unknownId, std::function<real_t(len_t,len_t,len_t)> *weightFunc)
        : EquationTerm(grid), unknownId(unknownId), weightFunc(weightFunc) { 
    weights = new real_t[grid->GetNCells()];
}

/**
 * Destructor
 */
WeightedTransientTerm::~WeightedTransientTerm(){
    delete [] weights;
}


/**
 * Rebuild the transient term.
 *
 * dt: Length of next time step to take.
 */
void WeightedTransientTerm::Rebuild(const real_t, const real_t dt, UnknownQuantityHandler *uqty) {
    this->dt = dt;
    this->xn = uqty->GetUnknownDataPrevious(this->unknownId);
    std::function<real_t(len_t,len_t,len_t)> func = *weightFunc;
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++)
        for(len_t i = 0; i < n1[ir]; i++)
            for(len_t j = 0; j < n2[ir]; j++){
                weights[offset + n1[ir]*j + i] = func(ir,i,j);
                offset += n1[ir]*n2[ir];
            }
 
}

/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 * x:       Value of the unknown quantity.
 */
void WeightedTransientTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *mat, const real_t* /*x*/
) {
    if (uqtyId == derivId)
        this->SetMatrixElements(mat, nullptr);
}

/**
 * Set the matrix elements corresponding to this
 * transient term.
 *
 * This term assumes that the linearized matrix equation
 * to solve is of the form
 *
 *   df/dt + Mf = -S
 *
 * where M is a linear matrix operator represented by 'mat',
 * and 'S' is a source term stored in 'rhs'.
 *
 * mat: Matrix to set elements of.
 * rhs: Equation RHS.
 */
void WeightedTransientTerm::SetMatrixElements(Matrix *mat, real_t *rhs) {
    const len_t N = grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, weights[i]/this->dt, ADD_VALUES);

    if (rhs != nullptr)
        for (len_t i = 0; i < N; i++)
            rhs[i] += weights[i]*this->xn[i] / this->dt;
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution (unused).
 */
void WeightedTransientTerm::SetVectorElements(real_t *vec, const real_t *xnp1) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] += weights[i]*(xnp1[i] - xn[i]) / this->dt;
}

