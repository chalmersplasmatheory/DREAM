/**
 * Implementation the EquationSystem class.
 */

#include <algorithm>
#include <vector>
#include <string>
#include <softlib/Timer.h>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IO.hpp"
#include "FVM/QuantityData.hpp"


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
EquationSystem::EquationSystem(
    FVM::Grid *rgrid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) : fluidGrid(rgrid), hottailGrid(hottailGrid), runawayGrid(runawayGrid) {
    
}

/**
 * Destructor.
 */
EquationSystem::~EquationSystem() {
    if (this->solver != nullptr)
        delete this->solver;
}

/**
 * Processes the system after it has been initialized and
 * prepares it for being solved.
 */
void EquationSystem::ProcessSystem() {
    // Construct a list of the unknowns that will appear
    // in any matrices later on
    len_t totsize = 0;
    const len_t N = unknowns.GetNUnknowns();
    for (len_t i = 0; i < N; i++) {
        if (!unknown_equations[i]->IsPrescribed()) {
            nontrivial_unknowns.push_back(i);
            totsize = unknowns[i]->GetGrid()->GetNCells();
        }
    }

    solver->Initialize(totsize, nontrivial_unknowns);
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
        //unknown_equations.insert(blockrow, new UnknownQuantityEquation(GetUnknown(blockrow)));
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
 * Solve this equation system.
 */
void EquationSystem::Solve() {
    this->currentTime = 0;

    // TODO Set initial state (or ensure that it has been set?)
    
    Timer tim;
    while (!timestepper->IsFinished(currentTime)) {
        real_t dt = timestepper->NextStep(currentTime);

        solver->Solve(currentTime, dt);
    }

    string duration = tim.ToString();

    DREAM::IO::PrintInfo("Solved equation system in %s.", duration.c_str());
}

