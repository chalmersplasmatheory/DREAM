/**
 * Implementation of an Euler backward transient term multiplied
 * by an arbitrary grid-dependent weight function. 
 * Evaluation of weights must be implemented in derived classes. 
 */

#include <iostream>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Equation/LinearTransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM::FVM;


/**
 * Constructor.
 */
LinearTransientTerm::LinearTransientTerm(Grid *grid, const len_t unknownId)
        : DiagonalTerm(grid), unknownId(unknownId){}


/**
 * Rebuild the transient term.
 *
 * dt: Length of next time step to take.
 */
void LinearTransientTerm::Rebuild(const real_t, const real_t dt, UnknownQuantityHandler *uqty) {
    this->dt = dt;
    this->xn = uqty->GetUnknownDataPrevious(this->unknownId);

    if(!hasBeenInitialized){
        InitializeWeights();
        hasBeenInitialized = true;
    }
    if(TermDependsOnUnknowns()) 
        SetWeights();
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
void LinearTransientTerm::SetMatrixElements(Matrix *mat, real_t *rhs) {
	// Transient term disabled? (mainly applicable for
	// equation system initializer, which seeks solutions for
	// which dX/dt -> 0)
	if (this->dt == 0)
		return;

    const len_t N = grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, weights[i]/this->dt, ADD_VALUES);

    if (rhs != nullptr)
        for (len_t i = 0; i < N; i++)
            rhs[i] -= weights[i]*this->xn[i] / this->dt;
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution (unused).
 */
void LinearTransientTerm::SetVectorElements(real_t *vec, const real_t *xnp1) {
	// Transient term disabled? (mainly applicable for
	// equation system initializer, which seeks solutions for
	// which dX/dt -> 0)
	if (this->dt == 0)
		return;

    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] += weights[i]*(xnp1[i] - xn[i]) / this->dt;
}

