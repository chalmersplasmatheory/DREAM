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
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/OutputGeneratorSFile.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
SolverLinearlyImplicit::SolverLinearlyImplicit(
    FVM::UnknownQuantityHandler *unknowns, 
    vector<UnknownQuantityEquation*> *unknown_equations,
    EquationSystem *eqsys,
    enum OptionConstants::linear_solver ls
) : Solver(unknowns, unknown_equations, eqsys, ls) {

    this->timeKeeper = new FVM::TimeKeeper("Solver linear");
    this->timerTot = this->timeKeeper->AddTimer("total", "Total time");
    this->timerRebuild = this->timeKeeper->AddTimer("rebuildtot", "Rebuild coefficients");
    this->timerMatrix = this->timeKeeper->AddTimer("matrix", "Construct matrix");
    this->timerInvert = this->timeKeeper->AddTimer("invert", "Invert matrix");
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

    // Select linear solver
    this->SelectLinearSolver(size);

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t id = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = this->unknown_equations->at(id);

        unknownToMatrixMapping[id] = 
            matrix->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros(), id);
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
void SolverLinearlyImplicit::SetInitialGuess(const real_t* /*guess*/) {
    /*if (guess != nullptr) {
        PetscScalar *x0;
        VecGetArray(petsc_sol, &x0);

        for (len_t i = 0; i < this->matrix_size; i++)
            x0[i] = guess[i];

        VecRestoreArray(petsc_sol, &x0);
    }*/
    // The initial guess is taken from the UnknownQuantityHandler,
    // and so this routine is not necessary...
}

/**
 * Solve the system of equations.
 *
 * t:  Time at which to solve the system.
 * dt: Time step to take.
 */
void SolverLinearlyImplicit::Solve(const real_t t, const real_t dt) {
    this->nTimeStep++;

    this->timeKeeper->StartTimer(timerTot);

    this->timeKeeper->StartTimer(timerRebuild);
    RebuildTerms(t, dt);
    this->timeKeeper->StopTimer(timerRebuild);

    real_t *S;
    VecGetArray(petsc_S, &S);
    this->timeKeeper->StartTimer(timerMatrix);
    BuildMatrix(t, dt, matrix, S);
    this->timeKeeper->StopTimer(timerMatrix);

    // Negate vector
    // We do this since in DREAM, we write the equation as
    //
    //   Mx + S = 0
    //
    // whereas PETSc solves the equation
    //
    //   Ax = b
    //
    // Thus, b = -S
    for (len_t i = 0; i < matrix->GetNRows(); i++)
        S[i] = -S[i];

    this->SaveDebugInfo(this->nTimeStep, matrix, S);

    VecRestoreArray(petsc_S, &S);

    // Apply preconditioner (if enabled)
    this->Precondition(matrix, petsc_S);

    this->timeKeeper->StartTimer(timerInvert);
    inverter->Invert(matrix, &petsc_S, &petsc_S);
    this->timeKeeper->StopTimer(timerInvert);

    // Undo preconditioner (if enabled)
    this->UnPrecondition(petsc_S);

    // Store solution
    unknowns->Store(this->nontrivial_unknowns, petsc_S);

    this->timeKeeper->StopTimer(timerTot);

    this->IterationFinished();
}

/**
 * Print timing information for this solver.
 */
void SolverLinearlyImplicit::PrintTimings() {
    this->timeKeeper->PrintTimings(true, 0);
    this->Solver::PrintTimings_rebuild();
}

/**
 * Save timing information to the given SFile object.
 *
 * sf:   SFile object to save timing information to.
 * path: Path in file to save timing information to.
 */
void SolverLinearlyImplicit::SaveTimings(SFile *sf, const string& path) {
    this->timeKeeper->SaveTimings(sf, path);

    sf->CreateStruct(path+"/rebuild");
    this->Solver::SaveTimings_rebuild(sf, path+"/rebuild");
}

/**
 * Save debug information, if enabled.
 *
 * it:  Time step index.
 * mat: Linear operator matrix.
 * rhs: Right-hand side vector.
 */
void SolverLinearlyImplicit::SaveDebugInfo(len_t it, FVM::Matrix *mat, const real_t *rhs) {
    if (this->savetimestep == it || this->savetimestep == 0) {
        string suffix = "_" + to_string(it);

        if (this->savematrix) {
            string matname;
            if (this->savetimestep == 0)
                matname = "petsc_mat" + suffix;
            else
                matname = "petsc_mat";

            mat->View(FVM::Matrix::BINARY_MATLAB, matname);
        }

        if (this->saverhs) {
            string rhsname;
            if (this->savetimestep == 0)
                rhsname = "rhs" + suffix + ".mat";
            else
                rhsname = "rhs.mat";

            SFile *sf = SFile::Create(rhsname, SFILE_MODE_WRITE);
            sf->WriteList("rhs", rhs, mat->GetNRows());
            sf->Close();
        }

        // Save full output?
        if (this->savesystem) {
            string outname = "debugout";
            if (this->savetimestep == 0)
                outname += suffix;
            outname += ".h5";

            OutputGeneratorSFile *outgen = new OutputGeneratorSFile(this->eqsys, outname, true);
            outgen->SaveCurrent();
            delete outgen;
        }

        if (this->printmatrixinfo)
            mat->PrintInfo();
    }
}

/**
 * Enable or disable debug mode (i.e. writing linear operator matrix
 * to file)
 *
 * printinfo:  If true, prints matrix debug info in every time step.
 * savematrix: If true, saves the linear operator matrix using a PETSc Viewer.
 * saverhs:    If true, saves the RHS vector to a MAT file.
 * timestep:   Index of time step for which to save the matrix. If '0', saves the
 *             matrix in all time steps.
 * iteration:  Index of iteration for which to save the matrix.
 * savesystem: If true, saves the full equation system, including grid information,
 *             to a proper DREAMOutput file. However, only the most recently obtained
 *             solution is saved.
 */
void SolverLinearlyImplicit::SetDebugMode(
    bool printinfo, bool savematrix, bool saverhs, int_t timestep,
    bool savesystem
) {
    this->printmatrixinfo = printinfo;
    this->savematrix = savematrix;
    this->saverhs = saverhs;
    this->savetimestep = timestep;
    this->savesystem = savesystem;
}

/**
 * Write basic data/statistics from this solver.
 *
 * sf:   SFile object to use for writing.
 * name: Name of group within file to store data in.
 */
void SolverLinearlyImplicit::WriteDataSFile(SFile *sf, const std::string &name) {
    sf->CreateStruct(name);

    int32_t type = (int32_t)OptionConstants::SOLVER_TYPE_LINEARLY_IMPLICIT;
    sf->WriteList(name+"/type", &type, 1);
}

