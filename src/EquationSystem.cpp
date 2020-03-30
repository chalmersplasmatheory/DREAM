/**
 * Implementation the EquationSystem class.
 */

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
 * Add an unknown quantity to the equation system.
 *
 * name: Name of unknown quantity.
 * grid: Grid on which the quantity is defined.
 */
int_t EquationSystem::SetUnknown(const string& name, FVM::Grid *grid) {
    struct unknown_qty *u = new struct unknown_qty(name, grid);
    unknowns.push_back(u);

    // Return ID of quantity
    return (unknowns.size()-1);
}

