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
    vector<UnknownQuantityEquation*> *unknown_equations,
    const PetscInt maxiter, const real_t reltol,
    bool verbose
) : Solver(unknowns, unknown_equations), maxiter(maxiter), reltol(reltol), verbose(verbose) {
    
    this->x_2norm  = new real_t[unknown_equations->size()];
    this->dx_2norm = new real_t[unknown_equations->size()];
}

/**
 * Destructor.
 */
SolverSNES::~SolverSNES() {
    if (jacobian != nullptr)
        delete jacobian;

    delete [] this->dx_2norm;
    delete [] this->x_2norm;

    VecDestroy(&this->petsc_F);
    VecDestroy(&this->petsc_sol);
    SNESDestroy(&this->snes);
}

/**
 * Calculate the 2-norm of the given vector separately for each
 * non-trivial unknown in the equation system. Thus, if there are
 * N non-trivial unknowns in the equation system, the output vector
 * 'retvec' will contain N elements, each holding the 2-norm of the
 * corresponding section of the input vector 'vec' (which is, for
 * example, a solution vector).
 */
void SolverSNES::CalculateNonTrivial2Norm(const real_t *vec, real_t *retvec) {
    len_t offset = 0, i = 0;
    for (auto id : this->nontrivial_unknowns) {
        FVM::UnknownQuantity *uqn = this->unknowns->GetUnknown(id);
        const len_t N = uqn->NumberOfElements();

        retvec[i] = 0;
        for (len_t j = 0; j < N; j++)
            retvec[i] += vec[offset+j]*vec[offset+j];

        retvec[i] = sqrt(retvec[i]);

        offset += N;
        i++;
    }
}

/**
 * Returns the name of the specified non-trivial unknown quantity.
 *
 * idx: Index into 'this->nontrivial_unknowns' of the non-trivial unknown
 *      to return the name of.
 */
const string& SolverSNES::GetNonTrivialName(const len_t idx) {
    return this->unknowns->GetUnknown(this->nontrivial_unknowns[idx])->GetName();
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
    // Set function evaluator
    SNESSetFunction(snes, petsc_F, &SNES_set_function, this);
    // Set jacobian evaluator
    SNESSetJacobian(snes, jacobian->mat(), jacobian->mat(), &SNES_set_jacobian, this);
    // Set solution monitor
    SNESMonitorSet(snes, &SNES_solution_obtained, this, nullptr);
    // Set convergence test
    SNESSetConvergenceTest(snes, &SNES_convergence_test, this, nullptr);

    // Newton's method with line search
    SNESSetType(snes, SNESNEWTONLS);
    SNESSetTolerances(snes,
        PETSC_DEFAULT,      // abstol
        PETSC_DEFAULT,      // reltol
        PETSC_DEFAULT,      // stol  (|| dx || < stol*|| x ||)
        this->MaxIter(),    // max iterations
        -1                  // max function evaluations
    );

    // Line-search options
    SNESLineSearch ls;
    SNESGetLineSearch(snes, &ls);

    SNESLineSearchSetType(ls, SNESLINESEARCHBASIC);
    //SNESLineSearchSetDamping(ls, .28);
    SNESLineSearchSetDamping(ls, .27);

    SNESLineSearchSetTolerances(ls,
        PETSC_DEFAULT,      // steptol
        1e50,               // maxstep
        PETSC_DEFAULT,      // reltol
        PETSC_DEFAULT,      // abstol
        PETSC_DEFAULT,      // lambda tolerance
        PETSC_DEFAULT       // max iterations
    );
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
    /*if (iteration == 8) {
        this->jacobian->View(FVM::Matrix::BINARY_MATLAB, "petsc_jacobian");
        this->_EvaluateJacobianNumerically(this->jacobian);

        throw SolverException("I want to stop now!");
    }*/
    //if (iteration == 1) {
        /*SFile *sf = SFile::Create("vectors/vector" + std::to_string(iteration) + ".mat", SFILE_MODE_WRITE);

        PetscInt size;
        VecGetSize(petsc_sol, &size);

        PetscScalar *x_arr = new PetscScalar[size];
        PetscScalar *dx_arr = new PetscScalar[size];
        PetscScalar *F_arr = new PetscScalar[size];
        PetscInt *idx = new PetscInt[size];
        for (PetscInt i = 0; i < size; i++)
            idx[i] = i;

        VecGetValues(petsc_F,   size, idx, F_arr);

        Vec x, dx;
        SNESGetSolution(this->snes, &x);
        SNESGetSolutionUpdate(this->snes, &dx);
        VecGetValues(x, size, idx, x_arr);
        VecGetValues(dx, size, idx, dx_arr);

        sf->WriteList("F", F_arr, size);
        sf->WriteList("dx", dx_arr, size);
        sf->WriteList("x", x_arr, size);
        sf->Close();

        delete [] F_arr;
        delete [] dx_arr;
        delete [] x_arr;
        delete [] idx;*/
    //}
}

/**
 * Solve the given stage for the non-linear
 * equation system.
 *
 * stage: Stage of the equation system to solve.
 *
 * RETURNS 'true' if the solve succeeded (i.e. did not
 * result in a diverging solution).
 */
void SolverSNES::Solve(const real_t t, const real_t dt) {
    this->t  = t;
    this->dt = dt;

    // Run SNES
    SNESSolve(this->snes, NULL, this->petsc_sol);

    SNESConvergedReason reason;
    SNESGetConvergedReason(this->snes, &reason);

    // Ensure that we're not diverging
    // (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html)
    if (reason < 0) {
        string rsn;

        switch (reason) {
            case SNES_DIVERGED_FUNCTION_DOMAIN: rsn = "the new x location passed to the function is not in the domain of F"; break;
            case SNES_DIVERGED_FUNCTION_COUNT:  rsn = "the user provided function has been called more times than the final argument to SNESSetTolerances()"; break;
            case SNES_DIVERGED_LINEAR_SOLVE:    rsn = "the linear solve failed"; break;
            case SNES_DIVERGED_FNORM_NAN:       rsn = "the norm of F is NaN"; break;
            case SNES_DIVERGED_MAX_IT:          rsn = "the maximum number of iterations reached without convergence"; break;
            case SNES_DIVERGED_LINE_SEARCH:     rsn = "the line search failed"; break;
            case SNES_DIVERGED_INNER:           rsn = "inner solve failed"; break;
            case SNES_DIVERGED_LOCAL_MIN:       rsn = "|| J^T b || is small, implies converged to local minimum of F()"; break;
            case SNES_DIVERGED_DTOL:            rsn = "|| F || > divtol*|| F_initial ||"; break;
            case SNES_DIVERGED_JACOBIAN_DOMAIN: rsn = "Jacobian calculation does not make sense"; break;
            case SNES_DIVERGED_TR_DELTA:        rsn = "TR_DELTA"; break;

            default:
                throw SolverException(
                    "The non-linear (SNES) solver failed to converge with reason: %d.",
                    reason
                );
        }

        throw SolverException(
            "The non-linear (SNES) solver failed to converge with reason: %s.",
            rsn.c_str()
        );
    } /*else
        throw SolverException("I want to quit now.");*/
}

