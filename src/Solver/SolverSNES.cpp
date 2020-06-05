/**
 * Interface to the Scalable Non-linear Equation Solver (SNES, part of PETSc).
 */

#include <vector>
#include <softlib/SFile.h>
#include "DREAM/Solver/SolverSNES.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
SolverSNES::SolverSNES(
    FVM::UnknownQuantityHandler *unknowns, 
    vector<UnknownQuantityEquation*> *unknown_equations
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
        len_t id = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = this->unknown_equations->at(id);

        unknownToMatrixMapping[id] =
            jacobian->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros_jac());
    }

    jacobian->ConstructSystem();

    // Construct solution and function vectors
    VecCreateSeq(PETSC_COMM_WORLD, size, &this->petsc_F);
    VecDuplicate(this->petsc_F, &this->petsc_sol);

    // Initialize SNES
    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFunction(snes, petsc_F, &SNES_set_function, this);
    SNESSetJacobian(snes, jacobian->mat(), jacobian->mat(), &SNES_set_jacobian, this);
    SNESMonitorSet(snes, &SNES_solution_obtained, this, nullptr);

    // Newton's method with line search
    SNESSetType(snes, SNESNEWTONLS);
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
void SolverSNES::StoreSolution(len_t iteration) {
    unknowns->Store(nontrivial_unknowns, petsc_sol);

    // DEBUG
    if (iteration == 1) {
        SFile *sf = SFile::Create("vector.mat", SFILE_MODE_WRITE);

        PetscInt size;
        VecGetSize(petsc_sol, &size);

        PetscScalar *x_arr = new PetscScalar[size];
        PetscScalar *F_arr = new PetscScalar[size];
        PetscInt *idx = new PetscInt[size];
        for (PetscInt i = 0; i < size; i++)
            idx[i] = i;

        VecGetValues(petsc_sol, size, idx, x_arr);
        VecGetValues(petsc_F,   size, idx, F_arr);

        sf->WriteList("F", F_arr, size);
        sf->WriteList("x", x_arr, size);
        sf->Close();

        delete [] F_arr;
        delete [] x_arr;
        delete [] idx;
    }
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

    // DEBUG
    /*PetscInt idx = 6;
    PetscScalar v;
    VecGetValues(petsc_sol, 1, &idx, &v);
    printf("f0 = %e\n", v);*/

    // Run SNES
    SNESSolve(this->snes, NULL, this->petsc_sol);
}

