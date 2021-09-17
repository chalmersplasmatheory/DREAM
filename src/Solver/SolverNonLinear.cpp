/**
 * Implementation of a custom Newton solver which only utilizes
 * the linear solvers of PETSc.
 */
#include <iostream>

#include "DREAM/IO.hpp"
#include "DREAM/OutputGeneratorSFile.hpp"
#include "DREAM/Solver/Backtracker.hpp"
#include "DREAM/Solver/PhysicalStepAdjuster.hpp"
#include "DREAM/Solver/SolverNonLinear.hpp"
#include <string>
#include <vector>


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
SolverNonLinear::SolverNonLinear(
    FVM::UnknownQuantityHandler *unknowns,
    vector<UnknownQuantityEquation*> *unknown_equations,
    EquationSystem *eqsys,
    enum OptionConstants::linear_solver ls,
    enum OptionConstants::linear_solver bk,
    enum OptionConstants::newton_step_adjuster nsa,
    const int_t maxiter, const real_t reltol,
    bool verbose
) : Solver(unknowns, unknown_equations, ls, bk), eqsys(eqsys),
    maxiter(maxiter), reltol(reltol), verbose(verbose), stepAdjusterType(nsa) {

    this->timeKeeper = new FVM::TimeKeeper("Solver non-linear");
    this->timerTot = this->timeKeeper->AddTimer("total", "Total time");
    this->timerRebuild = this->timeKeeper->AddTimer("rebuildtot", "Rebuild coefficients");
    this->timerResidual = this->timeKeeper->AddTimer("residual", "Construct residual");
    this->timerJacobian = this->timeKeeper->AddTimer("jacobian", "Construct jacobian");
    this->timerInvert = this->timeKeeper->AddTimer("invert", "Invert jacobian");
}

/**
 * Destructor.
 */
SolverNonLinear::~SolverNonLinear() {
    Deallocate();

    delete this->timeKeeper;
}


/**
 * "Accept" the current solution and prepare for taking
 * another Newton step.
 */
void SolverNonLinear::AcceptSolution() {
    real_t *x = this->x1;
    this->x1  = this->x0;
    this->x0  = x;

    this->StoreSolution(x);
}

/**
 * Allocate memory for all objects used by this solver.
 */
void SolverNonLinear::Allocate() {
    this->AllocateJacobianMatrix();

    const len_t N = jacobian->GetNRows();

    // Select linear solver
    this->SelectLinearSolver(N);

    VecCreateSeq(PETSC_COMM_WORLD, N, &this->petsc_F);
    VecCreateSeq(PETSC_COMM_WORLD, N, &this->petsc_dx);

    this->x0 = new real_t[N];
    this->x1 = new real_t[N];
    this->dx = new real_t[N];
    this->xinit = new real_t[N];

    this->x_2norm  = new real_t[this->unknown_equations->size()];
    this->dx_2norm = new real_t[this->unknown_equations->size()];
}

/**
 * Allocates memory for and properly sets up the jacobian matrix.
 * If the jacobian matrix has previously been allocated, it will
 * first be deleted.
 *
 * WHY DO WE CALL THIS METHOD MORE THAN ONCE?
 * In the very first iteration, many elements of the jacobian matrix are often
 * identically zero. If we don't insert the zeros explicitly, PETSc will remove
 * the memory we allocated for them on the first call to 'Assemble()' requiring
 * the memory to be reallocated in the next iteration (which may take a _very_
 * long time). However, if we insert the zeros explicitly, the linear solver
 * will not be able to tell that the elements are in fact non-zero and will
 * take ages to solve the system. As a compromise, we would like to use the
 * non-zero pattern obtained in the second Newton iteration of the simulation,
 * we should be very close to the non-zero pattern of the remainder of the
 * simulation. Hence, we call this method after the first iteration is finished
 * to completely reset the matrix, including the non-zero pattern. Reallocating
 * all the memory in one go is significantly faster than asking PETSc to
 * allocate memory specifically for all the elements we would like to add in the
 * second iteration.
 */
void SolverNonLinear::AllocateJacobianMatrix() {
    if (this->jacobian != nullptr)
        delete this->jacobian;
    this->jacobian = new FVM::BlockMatrix();

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t id = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = this->unknown_equations->at(id);

        unknownToMatrixMapping[id] =
            this->jacobian->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros_jac(), id);
    }

    this->jacobian->ConstructSystem();
}

/**
 * Deallocate memory used by this solver.
 */
void SolverNonLinear::Deallocate() {
    if (backupInverter != nullptr)
        delete backupInverter;

    delete mainInverter;
    delete jacobian;

    delete [] this->x_2norm;
    delete [] this->dx_2norm;

    delete [] this->x0;
    delete [] this->x1;
    delete [] this->dx;
    delete [] this->xinit;

    VecDestroy(&this->petsc_F);
    VecDestroy(&this->petsc_dx);
}

/**
 * Returns the name of the specified non-trivial unknown quantity.
 *
 * idx: Index into 'this->nontrivial_unknowns' of the non-trivial unknown
 *      to return the name of.
 */
const string& SolverNonLinear::GetNonTrivialName(const len_t idx) {
    return this->unknowns->GetUnknown(this->nontrivial_unknowns[idx])->GetName();
}

/**
 * Initialize the solver.
 */
void SolverNonLinear::initialize_internal(
    const len_t, vector<len_t>&
) {
    this->Allocate();

    if (this->convChecker == nullptr)
        this->SetConvergenceChecker(
            new ConvergenceChecker(unknowns, this->nontrivial_unknowns, this->reltol)
        );

    InitStepAdjuster(this->stepAdjusterType);
}

/**
 * Construct and initialize the Newton step adjuster algorithm
 * to use with this Newton solver.
 */
void SolverNonLinear::InitStepAdjuster(enum OptionConstants::newton_step_adjuster nsa) {
    switch (nsa) {
        case OptionConstants::NEWTON_STEP_ADJUSTER_BACKTRACK:
            this->adjuster = new Backtracker(this->nontrivial_unknowns, this->unknowns, this->ionHandler);
            break;

        default:
            this->adjuster = new PhysicalStepAdjuster(this->nontrivial_unknowns, this->unknowns, this->ionHandler, this->Verbose());
            break;
    }
}

/**
 * Check if the solver has converged.
 */
bool SolverNonLinear::IsConverged(const real_t *x, const real_t *dx) {
    if (this->GetIteration() >= this->MaxIter()){
        throw SolverException(
            "Non-linear solver reached the maximum number of allowed "
            "iterations: " LEN_T_PRINTF_FMT ".",
            this->MaxIter()
        );
    }
    
    // always print verbose for the last few iterations before reaching max
    const len_t numVerboseBeforeMax = 3; 
    bool printVerbose = this->Verbose() || (this->MaxIter() - this->GetIteration())<=numVerboseBeforeMax;
    if (printVerbose)
        DREAM::IO::PrintInfo("ITERATION %d", this->GetIteration());

    return convChecker->IsConverged(x, dx, printVerbose);
}

/**
 * Set the initial guess for the solver.
 *
 * guess: Vector containing values of initial guess.
 */
void SolverNonLinear::SetInitialGuess(const real_t *guess) {
    if (guess != nullptr) {
        for (len_t i = 0; i < this->matrix_size; i++)
            this->x0[i] = guess[i];
    } else {
        for (len_t i = 0; i < this->matrix_size; i++)
            this->x0[i] = 0;
    }
}

/**
 * Revert the solution to the initial guess.
 */
void SolverNonLinear::ResetSolution() {
    this->unknowns->GetLongVectorPrevious(this->nontrivial_unknowns, this->x0);
    this->StoreSolution(this->x0);
}

/**
 * Solve the equation system (advance the system in time
 * by one step).
 *
 * t:  Time at which the current solution is given.
 * dt: Time step to take.
 *
 * (the obtained solution will correspond to time t'=t+dt)
 */
void SolverNonLinear::Solve(const real_t t, const real_t dt) {
    // Return to main matrix inverter (in case backup inverter
    // was used to complete last time step)
    this->SwitchToMainInverter();

    this->nTimeStep++;

    this->t  = t;
    this->dt = dt;

    this->timeKeeper->StartTimer(timerTot);

    try {
        this->_InternalSolve();
    } catch (FVM::FVMException &ex) {
        // Retry with backup-solver (if allowed and not already used)
        if (this->backupInverter != nullptr && this->inverter != this->backupInverter) {
            if (this->Verbose()) {
                DREAM::IO::PrintInfo(
                    "Main inverter failed to converge. Switching to backup inverter."
                );
                DREAM::IO::PrintError(ex.what());
            }

            // Retry solve
            this->SwitchToBackupInverter();
            this->_InternalSolve();
        } else  // Rethrow exception
            throw ex;
    }

    // Save basic statistics for step
    this->nIterations.push_back(this->iteration);
    this->usedBackupInverter.push_back(this->inverter == this->backupInverter);

    this->timeKeeper->StopTimer(timerTot);
}

void SolverNonLinear::_InternalSolve() {
    // Take Newton steps
    len_t iter = 0;
    const real_t *x, *dx;
    do {
        iter++;
        this->SetIteration(iter);

REDO_ITER:
        dx = this->TakeNewtonStep();
        // Solution rejected (solver likely switched)
        if (dx == nullptr) {
            if (iter < this->MaxIter())
                goto REDO_ITER;
            else
                throw SolverException("Maximum number of iterations reached while dx=nullptr.");
        }

        x  = UpdateSolution(dx);
        
        AcceptSolution();
    } while (!IsConverged(x, dx));
}

/**
 * Debugging routine for saving both the "analytically" computed
 * Jacobian, as well as the Jacobian evaluated numerically using
 * finite differences, to file. When this method is called, the
 * 'jacobian' variable is assumed to contain the "analytical"
 * Jacobian matrix for the current time/iteration. This routine
 * will then save that matrix, compute the corresponding numerical
 * Jacobian, and save that.
 *
 * name: Base name to use for files.
 */
void SolverNonLinear::SaveNumericalJacobian(const std::string& name) {
    this->_EvaluateJacobianNumerically(this->jacobian);
    this->jacobian->View(FVM::Matrix::BINARY_MATLAB, name + "_num");
    abort();
}

void SolverNonLinear::SaveJacobian() {
    this->jacobian->View(FVM::Matrix::BINARY_MATLAB, "petsc_jacobian");
}
void SolverNonLinear::SaveJacobian(const std::string& name) {
    this->jacobian->View(FVM::Matrix::BINARY_MATLAB, name);
}

/**
 * Store the current solution to the UnknownQuantityHandler.
 */
void SolverNonLinear::StoreSolution(const real_t *x) {
    this->unknowns->Store(this->nontrivial_unknowns, x);
}

/**
 * Calculate the next Newton step to take.
 */
const real_t *SolverNonLinear::TakeNewtonStep() {
    this->timeKeeper->StartTimer(timerRebuild);
    this->RebuildTerms(this->t, this->dt);
    this->timeKeeper->StopTimer(timerRebuild);

    // Evaluate function vector
    this->timeKeeper->StartTimer(timerResidual);
    real_t *fvec;
    VecGetArray(this->petsc_F, &fvec);
    this->BuildVector(this->t, this->dt, fvec, this->jacobian);
    VecRestoreArray(this->petsc_F, &fvec);
    this->timeKeeper->StopTimer(timerResidual);
    
    // Reconstruct the jacobian matrix after taking the first
    // iteration.
    // (See the comment above 'AllocateJacobianMatrix()' for
    // details about why we do this...)
    if (this->nTimeStep == 1 && this->iteration == 2)
        this->AllocateJacobianMatrix();

    // Evaluate jacobian
    this->timeKeeper->StartTimer(timerJacobian);
    this->BuildJacobian(this->t, this->dt, this->jacobian);
    this->timeKeeper->StopTimer(timerJacobian);


    // Print/save debug info (if requested)
    this->SaveDebugInfo(this->nTimeStep, this->iteration);

    // Apply preconditioner (if enabled)
    this->Precondition(this->jacobian, this->petsc_F);

    // Solve J*dx = F
    this->timeKeeper->StartTimer(timerInvert);
    inverter->Invert(this->jacobian, &this->petsc_F, &this->petsc_dx);

    if (inverter->GetReturnCode() != 0) {
        if (this->Verbose())
            DREAM::IO::PrintInfo("Switching to backup inverter... " INT_T_PRINTF_FMT, inverter->GetReturnCode());

        this->SwitchToBackupInverter();

        return nullptr;
    }

    this->timeKeeper->StopTimer(timerInvert);

    // Undo preconditioner (if enabled)
    this->UnPrecondition(this->petsc_dx);

    // Copy dx
    VecGetArray(this->petsc_dx, &fvec);
    for (len_t i = 0; i < this->matrix_size; i++)
        this->dx[i] = fvec[i];
    VecRestoreArray(this->petsc_dx, &fvec);
    
    return this->dx;
}

/**
 * Update the current solution with the Newton step 'dx'.
 *
 * dx: Newton step to take.
 */
const real_t *SolverNonLinear::UpdateSolution(const real_t *dx) {
    //real_t dampingFactor = MaximalPhysicalStepLength(x0,dx,iteration,nontrivial_unknowns,unknowns,ionHandler,id_uqn);
    real_t dampingFactor = this->adjuster->Adjust(iteration, x0, dx, this->petsc_F, this->jacobian);
    
    for (len_t i = 0; i < this->matrix_size; i++)
        this->x1[i] = this->x0[i] - dampingFactor*dx[i];
    
    return this->x1;
}

/**
 * Print timing information after the solve.
 */
void SolverNonLinear::PrintTimings() {
    this->timeKeeper->PrintTimings(true, 0);
    this->Solver::PrintTimings_rebuild();
}

/**
 * Save timing information to the given SFile object.
 *
 * sf:   SFile object to save timing information to.
 * path: Path in file to save timing information to.
 */
void SolverNonLinear::SaveTimings(SFile *sf, const string& path) {
    this->timeKeeper->SaveTimings(sf, path);

    sf->CreateStruct(path+"/rebuild");
    this->Solver::SaveTimings_rebuild(sf, path+"/rebuild");
}

/**
 * Save debugging information for the current iteration.
 *
 * iTimeStep:  Current time step index.
 * iIteration: Current iteration index.
 */
void SolverNonLinear::SaveDebugInfo(
    len_t iTimeStep, len_t iIteration
) {
    if ((this->savetimestep == iTimeStep &&
        (this->saveiteration == iIteration || this->saveiteration == 0)) ||
         this->savetimestep == 0) {

        string suffix = "_" + to_string(iTimeStep) + "_" + to_string(iIteration);
        
        if (this->savejacobian) {
            string jacname;
            if (this->savetimestep == 0 || this->saveiteration == 0)
                jacname = "petsc_jac" + suffix;
            else
                jacname = "petsc_jac";

            SaveJacobian(jacname);
        }

        if (this->savevector) {
            string resname;
            if (this->savetimestep == 0 || this->saveiteration == 0)
                resname = "residual" + suffix + ".mat";
            else
                resname = "residual.mat";

            real_t *fvec;
            VecGetArray(this->petsc_F, &fvec);

            SFile *sf = SFile::Create(resname, SFILE_MODE_WRITE);
            sf->WriteList("F", fvec, this->jacobian->GetNRows());
            sf->Close();

            VecRestoreArray(this->petsc_F, &fvec);
        }

        if (this->savenumjac) {
            string jacname;
            if (this->savetimestep == 0 || this->saveiteration == 0)
                jacname = "petsc_jac" + suffix;
            else
                jacname = "petsc_jac";

            SaveNumericalJacobian(jacname);
        }

        if (this->savesolution) {
            string solname = "solution_dx";
            if (this->savetimestep == 0 || this->saveiteration == 0)
                solname += suffix;
            solname += ".mat";

            real_t *xvec;
            VecGetArray(this->petsc_dx, &xvec);

            SFile *sf = SFile::Create(solname, SFILE_MODE_WRITE);
            sf->WriteList("dx", xvec, this->jacobian->GetNRows());
            sf->Close();

            VecRestoreArray(this->petsc_dx, &xvec);
        }

        // Save full output?
        if (this->savesystem) {
            string outname = "debugout";
            if (this->savetimestep == 0 || this->saveiteration == 0)
                outname += suffix;
            outname += ".h5";

            OutputGeneratorSFile *outgen = new OutputGeneratorSFile(this->eqsys, outname);
            outgen->SaveCurrent();
            delete outgen;
        }

        if (this->printjacobianinfo)
            this->jacobian->PrintInfo();
    }
}

/**
 * Enable or disable debug mode.
 *
 * printjacobianinfo: If true, prints detailed debug information about the
 *                    PETSc matrix for the jacobian in every iteration.
 * savejacobian:      If true, saves the jacobian using a PETSc viewer in
 *                    the specified time step(s).
 * savesolution:      If true, saves the solution vector in the specified
 *                    time step(s).
 * savevector:        If true, saves the residual vector in the specified
 *                    time step(s).
 * savenumjac:        If true, calculates the jacobian numerically and saves
 *                    it using a PETSc viewer in the specified time step(s).
 * timestep:          Time step index to save debug info for. If 0, saves
 *                    the information in every iteration of every time step.
 * iteration:         Iteration of specified time step to save debug info for.
 * savesystem:        If true, saves the full equation system, including grid information,
 *                    to a proper DREAMOutput file. However, only the most recently obtained
 *                    solution is saved.
 */
void SolverNonLinear::SetDebugMode(
    bool printjacobianinfo, bool savesolution, bool savejacobian, bool savevector,
    bool savenumjac, int_t timestep, int_t iteration, bool savesystem
) {
    this->printjacobianinfo = printjacobianinfo;
    this->savejacobian      = savejacobian;
    this->savesolution      = savesolution;
    this->savevector        = savevector;
    this->savenumjac        = savenumjac;
    this->savetimestep      = timestep;
    this->saveiteration     = iteration;
    this->savesystem        = savesystem;
}


/**
 * Override switch to backup inverter.
 */
void SolverNonLinear::SwitchToBackupInverter() {
    // Switch inverter to use
    this->Solver::SwitchToBackupInverter();

    // Restore solution to initial guess for time step
    this->ResetSolution();
}

/**
 * Write basic data from the solver to the output file.
 * This data is mainly statistics about the solution.
 *
 * sf:   SFile object to use for writing.
 * name: Name of group within file to store data in.
 */
void SolverNonLinear::WriteDataSFile(SFile *sf, const std::string& name) {
    sf->CreateStruct(name);

    int32_t type = (int32_t)OptionConstants::SOLVER_TYPE_NONLINEAR;
    sf->WriteList(name+"/type", &type, 1);

    // Number of iterations per time step
    sf->WriteList(name+"/iterations", this->nIterations.data(), this->nIterations.size());

    // Whether or not backup inverter was used for a given time step
    len_t nubi = this->usedBackupInverter.size();
    int32_t *ubi = new int32_t[nubi];
    for (len_t i = 0; i < nubi; i++)
        ubi[i] = this->usedBackupInverter[i] ? 1 : 0;

    sf->WriteList(name+"/backupinverter", ubi, nubi);
    delete [] ubi;
}

