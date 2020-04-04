/**
 * Implementation the EquationSystem class.
 */

#include <algorithm>
#include <vector>
#include <string>
#include "DREAM/EquationSystem.hpp"


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
    if (this->matrix != nullptr)
        delete this->matrix;
}

/**
 * Returns the most recent data for the specified
 * unknown quantity.
 *
 * qty: ID of quantity to get data of.
 */
real_t *EquationSystem::GetUnknownData(const len_t qty) {
    return unknowns[qty]->data->Get();
}

/**
 * Returns the ID of the named unknown.
 *
 * name: Name of unknown quantity to get ID of.
 */
len_t EquationSystem::GetUnknownID(const std::string& name) {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++) {
        if ((*it)->name == name)
            return (it-unknowns.begin());
    }

    throw EquationSystemException(
        "No unknown quantity with name '%s' exists in the equation system."
    );
}

/**
 * Set the equation for the specified block matrix.
 * 
 * blockrow: Index of the unknown which a solution to the
 *           equation yields.
 * blockcol: Index of of the unknown which the equation is
 *           to be applied to.
 * eqn:      Equation to insert.
 */
void EquationSystem::SetEquation(len_t blockrow, len_t blockcol, FVM::Equation *eqn) {
    unknowns.at(blockrow)->equations[blockcol] = eqn;
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
 * Add an unknown quantity to the equation system.
 *
 * name: Name of unknown quantity.
 * grid: Grid on which the quantity is defined.
 */
len_t EquationSystem::SetUnknown(const string& name, FVM::Grid *grid) {
    struct unknown_qty *u = new struct unknown_qty(name, grid);
    unknowns.push_back(u);

    // Return ID of quantity
    return (unknowns.size()-1);
}

