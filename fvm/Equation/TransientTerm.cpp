/**
 * Implementation of a Euler backward transient term.
 *
 * df   f_{n+1} - f_n
 * -- ~ -------------
 * dt        dt
 *
 */

#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
TransientTerm::TransientTerm(Grid *grid) : EquationTerm(grid) { }

/**
 * Rebuild the transient term.
 *
 * dt: Length of next time step to take.
 */
void TransientTerm::Rebuild(const real_t dt) {
    this->dt = dt;
}

/**
 * Set the matrix elements corresponding to this term.
 * NOTE: This term should be applied last, as it modifies
 * all other elements of the matrix!
 *
 * This term assumes that the linearized matrix equation
 * to solve is of the form
 *
 *   df/dt = Mf + S
 *
 * where M is a linear matrix operator represented by 'mat'.
 * This term then discretizes the equation as
 *
 *   (f_{n+1} - f_n) / dt = Mf_{n+1} + S   <==>
 *
 *   (I - dt M) f_{n+1} = f_n + dt S
 *
 * ----------
 *  mat: Matrix to set elements of.
 *  rhs: Equation RHS.
 */
void TransientTerm::SetMatrixElements(Matrix *mat, real_t *rhs) {
    mat->IMinusDtA(this->dt);
    
    const len_t N = grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        rhs[i] *= this->dt;
}

