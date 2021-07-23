/**
 * Implementation of common routines for the 'Solver' routines.
 */

#include <iostream>

#include <vector>
#include "DREAM/IO.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "DREAM/EquationSystem.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"

// Linear solvers
#include "FVM/Solvers/MILU.hpp"
#ifdef PETSC_HAVE_MKL_PARDISO
#   include "FVM/Solvers/MIMKL.hpp"
#endif
#include "FVM/Solvers/MIMUMPS.hpp"
#include "FVM/Solvers/MISuperLU.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
Solver::Solver(
    FVM::UnknownQuantityHandler *unknowns,
    vector<UnknownQuantityEquation*> *unknown_equations,
    EquationSystem *eqsys,
    enum OptionConstants::linear_solver lsolve,
    enum OptionConstants::linear_solver bksolve
)
    : unknowns(unknowns), unknown_equations(unknown_equations), eqsys(eqsys), linearSolver(lsolve), backupSolver(bksolve) {

    this->solver_timeKeeper = new FVM::TimeKeeper("Solver rebuild");
    this->timerTot = this->solver_timeKeeper->AddTimer("total", "Total time");
    this->timerCqh = this->solver_timeKeeper->AddTimer("collisionhandler", "Rebuild coll. handler");
    this->timerREFluid = this->solver_timeKeeper->AddTimer("refluid", "Rebuild RunawayFluid");
    this->timerSPIHandler = this->solver_timeKeeper->AddTimer("spihandler", "Rebuild SPIHandler");
    this->timerRebuildTerms = this->solver_timeKeeper->AddTimer("equations", "Rebuild terms");
}

/**
 * Destructor.
 */
Solver::~Solver() {
    delete this->solver_timeKeeper;
    delete this->convChecker;

    if (this->diag_prec != nullptr)
        delete this->diag_prec;
}

/**
 * Build a jacobian matrix for the equation system.
 *
 * t:    Time to build the jacobian matrix for.
 * dt:   Length of time step to take.
 * mat:  Matrix to use for storing the jacobian.
 */
void Solver::BuildJacobian(const real_t, const real_t, FVM::BlockMatrix *jac) {
    // Reset jacobian matrix
    jac->Zero();

    // Iterate over (non-trivial) unknowns (i.e. those which appear
    // in the matrix system), corresponding to blocks in F and
    // rows in the Jacobian matrix.
    for (len_t uqnId : nontrivial_unknowns) {
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);
        map<len_t, len_t>& utmm = this->unknownToMatrixMapping;
        len_t matUqnId = utmm[uqnId];
        // Iterate over each equation term
        len_t operatorId = 0;
        for (auto it = eqn->GetOperators().begin(); it != eqn->GetOperators().end(); it++) {
            const real_t *x = unknowns->GetUnknownData(it->first);
        
            // "Differentiate with respect to the unknowns which
            // appear in the matrix"
            //   d (F_uqnId) / d x_derivId
            for (len_t derivId : nontrivial_unknowns) {
                len_t matDerivId = utmm[derivId];
                jac->SelectSubEquation(matUqnId, matDerivId);

                // - in the equation for                           x_uqnId
                // - differentiate the operator that is applied to x_it
                // - with respect to                               x_derivId
                it->second->SetJacobianBlock(it->first, derivId, jac, x);
            }

            operatorId++;
        }

        //printf("operatorId = " LEN_T_PRINTF_FMT "\n", operatorId);
    }
    jac->PartialAssemble();

    // Apply boundary conditions which overwrite elements
    for (len_t uqnId : nontrivial_unknowns) {
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);
        map<len_t, len_t>& utmm = this->unknownToMatrixMapping;
        len_t matUqnId = utmm[uqnId];

        // Iterate over each equation
        for (auto it = eqn->GetOperators().begin(); it != eqn->GetOperators().end(); it++) {
            const real_t *x = unknowns->GetUnknownData(it->first);

            // "Differentiate with respect to the unknowns which
            // appear in the matrix"
            //   d (eqn_uqnId) / d x_derivId
            for (len_t derivId : nontrivial_unknowns) {
                len_t matDerivId = utmm[derivId];
                jac->SelectSubEquation(matUqnId, matDerivId);
                // For logic, see comment in the for-loop above
                it->second->SetJacobianBlockBC(it->first, derivId, jac, x);
            }
        }
    }
    jac->Assemble();
}

/**
 * Build a linear operator matrix for the equation system.
 *
 * t:    Time to build the jacobian matrix for.
 * dt:   Length of time step to take.
 * mat:  Matrix to use for storing the jacobian.
 * rhs:  Right-hand-side in equation.
 */
void Solver::BuildMatrix(const real_t, const real_t, FVM::BlockMatrix *mat, real_t *S) {
    // Reset matrix and rhs
    mat->Zero();
    for (len_t i = 0; i < matrix_size; i++)
        S[i] = 0;

    // Build matrix
    for (len_t uqnId : nontrivial_unknowns) {
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);
        map<len_t, len_t>& utmm = this->unknownToMatrixMapping;
        len_t matUqnId = utmm[uqnId];
        for (auto it = eqn->GetOperators().begin(); it != eqn->GetOperators().end(); it++) {
            if (utmm.find(it->first) != utmm.end()) {
                mat->SelectSubEquation(matUqnId, utmm[it->first]);
                PetscInt vecoffs = mat->GetOffset(matUqnId);
                it->second->SetMatrixElements(mat, S + vecoffs);

            // The unknown to which this operator should be applied is a
            // "trivial" unknown quantity, meaning it does not appear in the
            // equation system matrix. We therefore build it as part of the
            // RHS vector.
            } else {
                PetscInt vecoffs = mat->GetOffset(matUqnId);
                const real_t *data = unknowns->GetUnknownData(it->first);
                it->second->SetVectorElements(S + vecoffs, data);
            }
        }
    }

    mat->Assemble();
}

/**
 * Build a function vector for the equation system.
 *
 * t:   Time to build the function vector for.
 * dt:  Length of time step to take.
 * vec: Vector to store evaluated equations in.
 * jac: Associated jacobian matrix.
 */
void Solver::BuildVector(const real_t, const real_t, real_t *vec, FVM::BlockMatrix *jac) {
    // Reset function vector
    for (len_t i = 0; i < matrix_size; i++)
        vec[i] = 0;

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t uqn_id = nontrivial_unknowns[i];
        len_t vecoffset = jac->GetOffset(i);
        unknown_equations->at(uqn_id)->SetVectorElements(vec+vecoffset, unknowns);
    }
}

/**
 * Calculate the 2-norm of the given vector separately for each
 * non-trivial unknown in the equation system. Thus, if there are
 * N non-trivial unknowns in the equation system, the output vector
 * 'retvec' will contain N elements, each holding the 2-norm of the
 * corresponding section of the input vector 'vec' (which is, for
 * example, a solution vector).
 *
 * vec:    Solution vector to calculate 2-norm of.
 * retvec: Vector which will contain result on return. The vector
 *         must have the same number of elements as there are
 *         non-trivial unknown quantities in the equation system.
 */
void Solver::CalculateNonTrivial2Norm(const real_t *vec, real_t *retvec) {
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
 * Initialize this solver.
 *
 * size:     Number of elements in full unknown vector.
 *           (==> jacobian is of size 'size-by-size').
 * unknowns: List of indices of unknowns to include in the
 *           function vectors/matrices.
 */
void Solver::Initialize(const len_t size, vector<len_t>& unknowns) {
    this->matrix_size = size;

    // Copy list of non-trivial unknowns (those which will
    // appear in the matrices that are built later on)
    nontrivial_unknowns = unknowns;

    this->initialize_internal(size, unknowns);
}

/**
 * This method is called whenever the solver finishes an
 * iteration. It in turn calls all registered callback
 * functions.
 */
void Solver::IterationFinished() {
    for (auto f : this->callbacks_iterationFinished)
        (*f)(this->eqsys->GetSimulation());
}

/**
 * Register a function to call whenever a solver iteration
 * has finished.
 */
void Solver::RegisterCallback_IterationFinished(
    iteration_finished_func_t f
) {
    this->callbacks_iterationFinished.push_back(f);
}

/**
 * Rebuild all equation terms in the equation system for
 * the specified time.
 *
 * t:  Time for which to rebuild the equation system.
 * dt: Length of time step to take next.
 */
void Solver::RebuildTerms(const real_t t, const real_t dt) {
    solver_timeKeeper->StartTimer(timerTot);

    this->ionHandler->Rebuild();
    // Rebuild ionHandler, collision handlers and RunawayFluid

    solver_timeKeeper->StartTimer(timerCqh);
    if (this->cqh_hottail != nullptr)
        this->cqh_hottail->Rebuild();
    if (this->cqh_runaway != nullptr)
        this->cqh_runaway->Rebuild();
    solver_timeKeeper->StopTimer(timerCqh);

    solver_timeKeeper->StartTimer(timerREFluid);
    this->REFluid -> Rebuild();
    solver_timeKeeper->StopTimer(timerREFluid);

    solver_timeKeeper->StartTimer(timerRebuildTerms);
    // Update prescribed quantities
    const len_t N = unknowns->Size();
    for (len_t i = 0; i < N; i++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(i);
        UnknownQuantityEquation *eqn = unknown_equations->at(i);

        if (eqn->IsPredetermined()) {
            eqn->RebuildEquations(t, dt, unknowns);
            FVM::PredeterminedParameter *pp = eqn->GetPredetermined();
            uqty->Store(pp->GetData(), 0, true);
        }
    }

    solver_timeKeeper->StartTimer(timerSPIHandler);
    if(this->SPI!=nullptr){
        this->SPI-> Rebuild(dt);
    }
    solver_timeKeeper->StopTimer(timerSPIHandler);

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t uqnId = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);

        for (auto it = eqn->GetOperators().begin(); it != eqn->GetOperators().end(); it++) {
            it->second->RebuildTerms(t, dt, unknowns);
        }
    }

    solver_timeKeeper->StopTimer(timerRebuildTerms);
    solver_timeKeeper->StopTimer(timerTot);
}

/**
 * Precondition the given matrix and RHS vector. This can improve
 * conditioning of the equation system to solve and should be called
 * on each solve (if enabled).
 */
void Solver::Precondition(FVM::Matrix *mat, Vec rhs) {
    if (this->diag_prec == nullptr)
        return;

    this->diag_prec->RescaleMatrix(mat);
    this->diag_prec->RescaleRHSVector(rhs);
}

/**
 * Transform the given solution, obtained from a preconditioned
 * equation system, so that it uses the same normalizations as
 * the rest of the code.
 */
void Solver::UnPrecondition(Vec x) {
    if (this->diag_prec == nullptr)
        return;

    this->diag_prec->UnscaleUnknownVector(x);
}

/**
 * Print timing information for the 'Rebuild' stage of the solver.
 * This stage looks the same for all solvers, and so is conveniently
 * defined here in the base class.
 */
void Solver::PrintTimings_rebuild() {
    this->solver_timeKeeper->PrintTimings(true, 0);
}

/**
 * Save timing information for the 'Rebuild' stage of the solver
 * to the given SFile object.
 */
void Solver::SaveTimings_rebuild(SFile *sf, const std::string& path) {
    this->solver_timeKeeper->SaveTimings(sf, path);
}

/**
 * Select the linear solver to use.
 *
 * N: Number of rows (or columns) in matrix to invert.
 */
void Solver::SelectLinearSolver(const len_t N) {
    if (this->linearSolver == this->backupSolver)
        throw SolverException(
            "The main and backup linear solvers may not be the same."
        );

    this->mainInverter = this->ConstructLinearSolver(N, this->linearSolver);
    this->inverter = this->mainInverter;

    if (this->backupSolver != OptionConstants::LINEAR_SOLVER_NONE)
        this->backupInverter = this->ConstructLinearSolver(N, this->backupSolver);
}

FVM::MatrixInverter *Solver::ConstructLinearSolver(const len_t N, enum OptionConstants::linear_solver ls) {
    if (ls == OptionConstants::LINEAR_SOLVER_LU)
        return new FVM::MILU(N);
    else if (ls == OptionConstants::LINEAR_SOLVER_MKL) {
#ifdef PETSC_HAVE_MKL_PARDISO
        return new FVM::MIMKL(N);
#else
        throw SolverException(
            "Your version of PETSc does not include support for Intel MKL PARDISO. "
            "To use this linear solver you must recompile PETSc."
        );
#endif
    } else if (ls == OptionConstants::LINEAR_SOLVER_MUMPS) {
#ifdef PETSC_HAVE_MUMPS
        return new FVM::MIMUMPS(N);
#else
        throw SolverException(
            "Your version of PETSc does not include support for MUMPS. "
            "To use this linear solver you must recompile PETSc."
        );
#endif
    } else if (ls == OptionConstants::LINEAR_SOLVER_SUPERLU) {
#ifdef PETSC_HAVE_SUPERLU
        return new FVM::MISuperLU(N);
#else
        throw SolverException(
            "Your version of PETSc does not include support for SuperLU. "
            "To use this linear solver you must recompile PETSc."
        );
#endif
    } else
        throw SolverException(
            "Unrecognized linear solver specified: %d.", ls
        );
}

/**
 * Set the convergence checker to use for the linear solver.
 */
void Solver::SetConvergenceChecker(ConvergenceChecker *cc) {
    if (this->convChecker != nullptr)
        delete this->convChecker;

    this->convChecker = cc;
}

/**
 * Set the preconditioner to use for solver. If 'nullptr', no
 * preconditioning is done.
 */
void Solver::SetPreconditioner(DiagonalPreconditioner *dp) {
    if (this->diag_prec != nullptr)
        delete this->diag_prec;

    this->diag_prec = dp;
}

/**
 * Switch to using the backup inverter instead of the main inverter.
 */
void Solver::SwitchToBackupInverter() {
    if (this->inverter == this->backupInverter)
        throw DREAMException("Backup matrix inverter failed with PETSc error: " INT_T_PRINTF_FMT, inverter->GetReturnCode());

    this->inverter = this->backupInverter;
}

/**
 * Switch to using the main inverter.
 */
void Solver::SwitchToMainInverter() {
    this->inverter = this->mainInverter;
}

/**
 * Empty routine for writing solver data to output file.
 * This method should be overridden where needed.
 */
void Solver::WriteDataSFile(SFile*, const std::string&) {}

