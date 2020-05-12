/**
 * Implementation of the linearly implicit solution method of the non-linear
 * equation system. With this solver, we first assume that the non-linear
 * function F(x) can be decomposed into
 *
 *   F(x) ~ M(t,x) x(t) + S(t,x),
 *
 * where M(t,x) is an operator matrix, S(t,x) is a vector and x(t) is the
 * unknown vector. Next, we assume that M(t,x) and S(t,x) vary slowly with time
 * such that
 *
 *   M(t_{n+1}, x_{n+1}) ~ M(t_n, x_n),
 *   S(t_{n+1}, x_{n+1}) ~ S(t_n, x_n),
 *
 * where x_{n+1} = x(t_{n+1}). With this approximation, we may write the
 * originally non-linear equation system "F(x) = 0" as the system of linear
 * equations
 *   
 *   F(x_{n+1}) ~ M(x_n) x_{n+1} + S(x_n) = 0,
 *
 * which is solved by
 *
 *   x_{n+1} = -M(x_n)^{-1} S(x_n).
 *
 * Given the initial state 'x_n' of the system, we may thus straightforwardly
 * advance the system in time.
 */

#include <vector>
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"
#include "FVM/Solvers/MILU.hpp"


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
SolverLinearlyImplicit::SolverLinearlyImplicit(
    FVM::UnknownQuantityHandler *unknowns, 
    vector<UnknownQuantityEquation*> *unknown_equations
) : Solver(unknowns, unknown_equations) {
}

/**
 * Destructor.
 */
SolverLinearlyImplicit::~SolverLinearlyImplicit() {
    delete this->matrix;
    delete this->inverter;

    VecDestroy(&this->petsc_sol);
    VecDestroy(&this->petsc_S);
}

/**
 * Initialize the solver.
 */
void SolverLinearlyImplicit::initialize_internal(
    const len_t size, std::vector<len_t>&
) {
    this->matrix = new FVM::BlockMatrix();
    this->inverter = new FVM::MILU(size);

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t id = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = this->unknown_equations->at(id);

        unknownToMatrixMapping[id] = 
            matrix->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros());
    }

    matrix->ConstructSystem();

    VecCreateSeq(PETSC_COMM_WORLD, size, &this->petsc_S);
    VecCreateSeq(PETSC_COMM_WORLD, size, &this->petsc_sol);
}

/**
 * Set the initial guess for the linear solver.
 *
 * guess: Initial guess. If 'nullptr', uses the previous
 *        solution as the initial guess.
 */
void SolverLinearlyImplicit::SetInitialGuess(const real_t *guess) {
    if (guess != nullptr) {
        PetscScalar *x0;
        VecGetArray(petsc_sol, &x0);

        for (len_t i = 0; i < this->matrix_size; i++)
            x0[i] = guess[i];

        VecRestoreArray(petsc_sol, &x0);
    }
}

/**
 * Solve the system of equations.
 *
 * t:  Time at which to solve the system.
 * dt: Time step to take.
 */
void SolverLinearlyImplicit::Solve(const real_t t, const real_t dt) {
    RebuildTerms(t, dt);

    real_t *S;
    VecGetArray(petsc_S, &S);
    BuildMatrix(t, dt, matrix, S);
    if (t == 0) {
        SFile *sf = SFile::Create("vec.mat", SFILE_MODE_WRITE);
        sf->WriteList("vec", S, matrix->GetNRows());
        sf->Close();
    }
    VecRestoreArray(petsc_S, &S);

    //matrix->PrintInfo();
    if (t == 0)
        matrix->View(FVM::Matrix::BINARY_MATLAB);
    //matrix->View(FVM::Matrix::ASCII_MATLAB);
    inverter->Invert(matrix, &petsc_S, &petsc_S);

    // Store solution
    unknowns->Store(this->nontrivial_unknowns, petsc_S);
}

