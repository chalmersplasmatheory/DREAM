/**
 * Implementation of matrix invertor based on LU
 * factorization (i.e. direct inversion).
 */

#include <algorithm>
#include <functional>
#include <utility>
#include <petscmat.h>
#include <petscvec.h>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Solvers/MIMKL.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * n: Number of elements in solution vector.
 */
MIMKL::MIMKL(const len_t n, bool verbose, std::vector<PetscInt> schurIndices)
    : verbose(verbose), schurIndices(std::move(schurIndices)) {
    this->ksp = nullptr;
    if (this->schurIndices.empty())
        KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    this->xn = n;
}

/**
 * Destructor.
 */
MIMKL::~MIMKL() {
    KSPDestroy(&this->ksp);
}

#ifdef PETSC_HAVE_MKL_PARDISO
/**
 * Factor a fresh matrix so PARDISO always sees the index set before symbolic
 * analysis, including when the Jacobian's sparsity pattern changes.
 */
static PetscErrorCode SolveSchur(Mat factor, Mat A, Vec b, Vec x,
    const std::vector<PetscInt>& indices, bool verbose
) {
    const MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(A));
    PetscInt rows;
    IS schur = nullptr;
    MatFactorError factorError;

    PetscFunctionBeginUser;
    PetscCall(MatGetSize(A, &rows, nullptr));
    PetscCheck(indices.size() < static_cast<size_t>(rows) && indices.front() >= 0 &&
        indices.back() < rows &&
        std::adjacent_find(indices.begin(), indices.end(), std::greater_equal<PetscInt>()) == indices.end(),
        comm, PETSC_ERR_ARG_OUTOFRANGE, "Expected sorted, unique Schur indices in a proper subset of the rows");

    // Use the same one-thread setting as the ordinary DREAM MKL solver.
    PetscCall(MatMkl_PardisoSetCntl(factor, 65, 1));
    if (verbose)
        PetscCall(MatMkl_PardisoSetCntl(factor, 68, 1));

    PetscCall(ISCreateGeneral(comm, static_cast<PetscInt>(indices.size()), indices.data(), PETSC_COPY_VALUES, &schur));
    const PetscErrorCode setError = MatFactorSetSchurIS(factor, schur);
    const PetscErrorCode destroyError = ISDestroy(&schur);
    PetscCall(setError);
    PetscCall(destroyError);

    PetscCall(MatLUFactorSymbolic(factor, A, nullptr, nullptr, nullptr));
    PetscCall(MatLUFactorNumeric(factor, A, nullptr));
    PetscCall(MatFactorGetError(factor, &factorError));
    PetscCheck(factorError == MAT_FACTOR_NOERROR, comm, PETSC_ERR_NOT_CONVERGED,
        "PARDISO Schur factorization failed with reason %d", static_cast<int>(factorError));

    // PETSc solves the dense Schur system and expands the full correction.
    PetscCall(MatSolve(factor, b, x));
    PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

/**
 * Solves the linear equation system represented by
 *
 *   Ax = b
 *
 * where A is a matrix, and b and x are vectors. A
 * pointer is returned to the solution, x.
 *
 * A: Matrix of size m-by-n representing the linear system.
 * b: Right-hand-side vector containing n elements.
 * x: Solution vector. Contains solution on return. Must be
 *    of size n at least.
 */
#ifdef PETSC_HAVE_MKL_PARDISO
void MIMKL::Invert(Matrix *A, Vec *b, Vec *x) {
    if (!this->schurIndices.empty()) {
        Mat factor = nullptr;
        this->errorcode = MatGetFactor(A->mat(), MATSOLVERMKL_PARDISO, MAT_FACTOR_LU, &factor);
        if (this->errorcode == 0)
            this->errorcode = SolveSchur(factor, A->mat(), *b, *x, this->schurIndices, this->verbose);
        const PetscErrorCode destroyError = MatDestroy(&factor);
        if (this->errorcode == 0)
            this->errorcode = destroyError;
        return;
    }

    Mat F;
    PC pc;

    KSPSetOperators(this->ksp, A->mat(), A->mat());
    
    // Set direct LU factorization
    KSPGetPC(this->ksp, &pc);
    PCSetType(pc, PCLU);
    PCFactorSetMatSolverType(pc, MATSOLVERMKL_PARDISO);
    PCFactorSetUpMatSolverType(pc);
    PCFactorGetMatrix(pc, &F);
    KSPSetType(this->ksp, KSPPREONLY);

    // Force PARDISO to only use one CPU thread
    // (otherwise it will crash with an "insufficient memory"
    // error)
    MatMkl_PardisoSetCntl(F, 65, 1);

    // Maximum iterations for refinement
    MatMkl_PardisoSetCntl(F, 8, 10);

    if (this->verbose)
        MatMkl_PardisoSetCntl(F, 68, 1);

    // Solve
    this->errorcode = KSPSolve(this->ksp, *b, *x);

    if (this->errorcode != 0)
        PCPostSolve(pc, this->ksp);
#else
void MIMKL::Invert(Matrix*, Vec*, Vec*) {
#endif
}

