/**
 * Implementation the EquationSystem class.
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <softlib/Timer.h>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/QuantityData.hpp"

// DEBUG
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
EquationSystem::EquationSystem(
    FVM::Grid *emptygrid, FVM::Grid *rgrid,
    enum OptionConstants::momentumgrid_type ht_type, FVM::Grid *hottailGrid,
    enum OptionConstants::momentumgrid_type re_type, FVM::Grid *runawayGrid
) : scalarGrid(emptygrid), fluidGrid(rgrid),
    hottailGrid(hottailGrid), runawayGrid(runawayGrid),
    hottailGrid_type(ht_type), runawayGrid_type(re_type) {
    
    this->initializer = new EqsysInitializer(&this->unknowns, &this->unknown_equations);
}

/**
 * Destructor.
 */
EquationSystem::~EquationSystem() {
    if (this->ionHandler != nullptr)
        delete this->ionHandler;
    if (this->solver != nullptr)
        delete this->solver;
    if (this->timestepper != nullptr)
        delete this->timestepper;

    if (this->cqh_hottail != nullptr)
        delete this->cqh_hottail;
    if (this->cqh_runaway != nullptr)
        delete this->cqh_runaway;

    if (this->REFluid != nullptr)
        delete this->REFluid;
    if (this->postProcessor != nullptr)
        delete this->postProcessor;
}

/**
 * Processes the system after it has been initialized and
 * prepares it for being solved. This includes setting
 * initial values for those unknown quantities which do
 * not yet have initial values.
 *
 * t0: Time at which to set initial values.
 */
void EquationSystem::ProcessSystem(const real_t t0) {
    // Construct a list of the unknowns that will appear
    // in any matrices later on
    len_t totsize = 0;
    const len_t N = unknowns.GetNUnknowns();
    bool unknownMissing = false;
    for (len_t i = 0; i < N; i++) {
        if (unknown_equations[i] == nullptr) {
            DREAM::IO::PrintError("No equation has been declared for unknown '%s'", unknowns.GetUnknown(i)->GetName().c_str());
            unknownMissing = true;
        } else {
            if (!unknown_equations[i]->IsPredetermined()) {
                nontrivial_unknowns.push_back(i);
                totsize += unknowns[i]->NumberOfElements();
            }
        }
    }

    // Set initial values
    this->initializer->Execute(t0);

    if (unknownMissing)
        throw EquationSystemException("While processing equation system: Equations not declared for some unknowns.");

    this->matrix_size = totsize;
}

/**
 * Set one equation of the specified unknown.
 */
void EquationSystem::SetEquation(const len_t blockrow, const len_t blockcol, FVM::Equation *eqn, const std::string& desc) {
    // Verify that the list is sufficiently large
    if (unknown_equations.size() < blockrow+1)
        unknown_equations.resize(unknowns.Size(), nullptr);

    // Does the unknown have any equations yet? If not, create
    // first the equation container
    if (unknown_equations[blockrow] == nullptr)
        unknown_equations[blockrow] = new UnknownQuantityEquation(blockrow, GetUnknown(blockrow), desc);

    unknown_equations[blockrow]->SetEquation(blockcol, eqn);

    if (desc != "")
        unknown_equations[blockrow]->SetDescription(desc);
}

/**
 * Same as 'SetEquation(len_t, len_t, Equation*)', but specifies
 * the unknowns by name rather than by index.
 */
void EquationSystem::SetEquation(len_t blockrow, const std::string& blockcol, FVM::Equation *eqn, const std::string& desc) {
    SetEquation(blockrow, GetUnknownID(blockcol), eqn, desc);
}
void EquationSystem::SetEquation(const std::string& blockrow, len_t blockcol, FVM::Equation *eqn, const std::string& desc) {
    SetEquation(GetUnknownID(blockrow), blockcol, eqn, desc);
}
void EquationSystem::SetEquation(const std::string& blockrow, const std::string& blockcol, FVM::Equation *eqn, const std::string& desc) {
    SetEquation(GetUnknownID(blockrow), GetUnknownID(blockcol), eqn, desc);
}

/**
 * Set the initial value of the specified unknown quantity. If
 * the initial value has previously been specified, it is overwritten.
 *
 * id:  ID of unknown quantity.
 * val: Initial value of the quantity.
 * t0:  Initial time.
 */
void EquationSystem::SetInitialValue(const std::string& name, const real_t *val, const real_t t0) {
    this->SetInitialValue(this->unknowns.GetUnknownID(name), val, t0);
}
void EquationSystem::SetInitialValue(const len_t id, const real_t *val, const real_t t0) {
    this->unknowns.SetInitialValue(id, val, t0);
}

/**
 * Set the equation solver for this system.
 * Note that this method must be called only AFTER 'ProcessSystem()'
 * has been called.
 *
 * solver: Solver to assign
 */
void EquationSystem::SetSolver(Solver *solver) {
    if (this->matrix_size == 0)
        throw EquationSystemException("It appears that either 'EquationSystem::ProcessSystem()' has not been called prior to 'SetSolver()', or the matrix system is trivial (matrix_size = 0).");

    this->solver = solver;
    this->solver->Initialize(this->matrix_size, this->nontrivial_unknowns);
}

/**
 * Solve this equation system.
 */
void EquationSystem::Solve() {
    this->currentTime = 0;
    this->times.push_back(this->currentTime);

    this->PrintNonTrivialUnknowns();
    this->PrintTrivialUnknowns();

    // TODO Set initial state (or ensure that it has been set?)

    // Set initial guess in solver
    const real_t *guess = unknowns.GetLongVector(this->nontrivial_unknowns);
    solver->SetInitialGuess(guess);
    delete [] guess;
    
    cout << "Beginning time advance..." << endl;

    Timer tim;
    len_t istep = 0;
    while (!timestepper->IsFinished()) {
        // XXX Should this be done this early?
        // Rebuild collision quantity handlers
        //this->Rebuild();

        // Take step
        real_t dt = timestepper->NextTime() - this->currentTime;

        solver->Solve(this->currentTime, dt);
        this->currentTime = timestepper->CurrentTime();
        istep++;

        // DEBUG
        /*if (istep == 1) {
            const len_t Nsize = unknowns.GetLongVectorSize(this->nontrivial_unknowns);
            real_t *vec = new real_t[Nsize];
            FVM::BlockMatrix *mat = static_cast<SolverLinearlyImplicit*>(solver)->GetMatrix();
            solver->BuildVector(this->currentTime, dt, vec, mat);

            const real_t *x = unknowns.GetLongVector(this->nontrivial_unknowns);

            SFile *sf = SFile::Create("solution.mat", SFILE_MODE_WRITE);
            sf->WriteList("F", vec, Nsize);
            sf->WriteList("x", x, Nsize);
            sf->Close();

            delete [] vec;
            delete [] x;
        }*/

        // Post-process solution (should be done before saving any
        // time step)
        this->postProcessor->Process(this->currentTime);

        // true = Really save the step (if it's false, we just
        // indicate that we have taken another timestep). This
        // should only be true for time steps which we want to
        // push to the output file.
        unknowns.SaveStep(this->currentTime, true);
        this->times.push_back(this->currentTime);

        otherQuantityHandler->StoreAll(this->currentTime);

        cout << istep << "... ";
        if (istep % 10 == 0) cout << endl;
    }

    cout << endl;

    string duration = tim.ToString();

    DREAM::IO::PrintInfo("Solved equation system in %s.", duration.c_str());
}

/**
 * Rebuild quantities that need to be updated between iterations
 */
void EquationSystem::Rebuild(){
    if (this->cqh_hottail != nullptr)
        this->cqh_hottail->Rebuild();
    if (this->cqh_runaway != nullptr)
        this->cqh_runaway->Rebuild();
    bool useApproximateEceffMethod = false; // true if we want to use the approximate Eceff method instead, which is significantly faster but with ~10-20% inaccuracy
    this->REFluid -> Rebuild(useApproximateEceffMethod);

}

