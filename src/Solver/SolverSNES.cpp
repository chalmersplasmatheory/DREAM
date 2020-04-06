/**
 * Interface to the Scalable Non-linear Equation Solver (SNES, part of PETSc).
 */

#include <vector>
#include "DREAM/Solver/SolverSNES.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
SolverSNES::SolverSNES(
    FVM::UnknownQuantityHandler *unknowns, 
    vector<UnknownQuantityEquation*> *unknonw_equations
) : Solver(unknowns, unknown_equations) {
}

/**
 * Destructor.
 */
SolverSNES::~SolverSNES() {
    if (jacobian != nullptr)
        delete jacobian;

    VecDestroy(&this->petsc_F);
    VecDestroy(&this->petsc_sol);
    SNESDestroy(&this->snes);
}

/**
 * Initialize this solver.
 *
 * size:                Number of elements in full unknown vector.
 *                      (==> jacobian is of size 'size-by-size').
 * nontrivial_unknowns: List of indices of unknowns to include in the
 *                      function vectors/matrices.
 */
void SolverSNES::initialize_internal(const len_t size, vector<len_t>& nontrivial_unknowns) {
    jacobian = new FVM::BlockMatrix();

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        UnknownQuantityEquation *eqn = this->unknown_equations->at(i);

        jacobian->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros_jac());
    }

    jacobian->ConstructSystem();

    // Initialize SNES
    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFunction(snes, petsc_F, &SNES_set_function, this);
    SNESSetJacobian(snes, jacobian->mat(), jacobian->mat(), &SNES_set_jacobian, this);
    SNESMonitorSet(snes, &SNES_solution_obtained, this, nullptr);

    // Newton's method with line search
    SNESSetType(snes, SNESNEWTONLS);

    // Construct solution and function vectors
    VecCreateSeq(PETSC_COMM_WORLD, size, &this->petsc_F);
    VecDuplicate(this->petsc_F, &this->petsc_sol);
}

/**
 * Set the initial guess for the Newton solution.
 *
 * guess: Initial guess. If 'nullptr', uses the previous
 *        solution as the initial guess.
 */
void SolverSNES::SetInitialGuess(const real_t *guess) {
    if (guess != nullptr) {
        PetscScalar *x0;
        VecGetArray(petsc_sol, &x0);

        for (len_t i = 0; i < this->matrix_size; i++)
            x0[i] = guess[i];

        VecRestoreArray(petsc_sol, &x0);
    }
}

/**
 * Store the current solution to the UnknownQuantityHandler.
 */
void SolverSNES::StoreSolution() {
    unknowns->Store(nontrivial_unknowns, petsc_sol);
}

/**
 * Solve the given stage for the non-linear
 * equation system.
 *
 * stage: Stage of the equation system to solve.
 */
void SolverSNES::Solve(const real_t t, const real_t dt) {
    this->t  = t;
    this->dt = dt;

    // Run SNES
    SNESSolve(this->snes, NULL, this->petsc_sol);
}

