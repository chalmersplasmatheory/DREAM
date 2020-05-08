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


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
EquationSystem::EquationSystem(
    FVM::Grid *rgrid,
    enum OptionConstants::momentumgrid_type ht_type, FVM::Grid *hottailGrid,
    enum OptionConstants::momentumgrid_type re_type, FVM::Grid *runawayGrid
) : fluidGrid(rgrid),
    hottailGrid(hottailGrid), runawayGrid(runawayGrid),
    hottailGrid_type(ht_type), runawayGrid_type(re_type) {
    
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
}

/**
 * Processes the system after it has been initialized and
 * prepares it for being solved. This includes setting
 * initial values for those unknown quantities which do
 * not yet have initial values.
 */
void EquationSystem::ProcessSystem() {
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

            // Set initial value if not already set
            if (!unknowns[i]->HasInitialValue()) {
                const real_t t0 = 0, dt = 0;
                DREAM::IO::PrintInfo("Automatically setting initial value for '%s'...", unknowns[i]->GetName().c_str());

                // Handle evaluables automatically
                if (unknown_equations[i]->IsEvaluable()) {
                    const len_t nu = unknowns[i]->NumberOfElements();
                    real_t *vec = new real_t[nu];
                    for (len_t i = 0; i < nu; i++)
                        vec[i] = 0.0;

                    unknown_equations[i]->RebuildEquations(t0, dt, &unknowns);
                    unknown_equations[i]->Evaluate(vec, &unknowns);

                    SetInitialValue(i, vec, t0);

                    delete [] vec;
                } else
                    throw EquationSystemException(
                        "Unable to automatically set initial value for "
                        "non-predetermined quantity: '%s'...",
                        unknowns[i]->GetName().c_str()
                    );
            }
        }
    }

    if (unknownMissing)
        throw EquationSystemException("While processing equation system: Equations not declared for some unknowns.");

    this->matrix_size = totsize;
}

/**
 * Set one equation of the specified unknown.
 */
void EquationSystem::SetEquation(const len_t blockrow, const len_t blockcol, FVM::Equation *eqn) {
    // Verify that the list is sufficiently large
    if (unknown_equations.capacity() < blockrow+1)
        unknown_equations.resize(unknowns.Size(), nullptr);

    // Does the unknown have any equations yet? If not, create
    // first the equation container
    if (unknown_equations[blockrow] == nullptr)
        unknown_equations[blockrow] = new UnknownQuantityEquation(GetUnknown(blockrow));

    unknown_equations[blockrow]->SetEquation(blockcol, eqn);
}

/**
 * Same as 'SetEquation(len_t, len_t, Equation*)', but specifies
 * the unknowns by name rather than by index.
 */
void EquationSystem::SetEquation(len_t blockrow, const std::string& blockcol, FVM::Equation *eqn) {
    SetEquation(blockrow, GetUnknownID(blockcol), eqn);
}
void EquationSystem::SetEquation(const std::string& blockrow, len_t blockcol, FVM::Equation *eqn) {
    SetEquation(GetUnknownID(blockrow), blockcol, eqn);
}
void EquationSystem::SetEquation(const std::string& blockrow, const std::string& blockcol, FVM::Equation *eqn) {
    SetEquation(GetUnknownID(blockrow), GetUnknownID(blockcol), eqn);
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
    cout << "Beginning time advance..." << endl;

    Timer tim;
    len_t istep = 0;
    while (!timestepper->IsFinished(currentTime)) {
        this->Rebuild();
        real_t dt = timestepper->NextStep(currentTime);

        solver->Solve(currentTime, dt);
        this->currentTime += dt;
        istep++;

        unknowns.SaveStep(this->currentTime);
        this->times.push_back(this->currentTime);

        cout << istep << "... ";
        if (istep % 10 == 0) cout << endl;
    }

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
}

