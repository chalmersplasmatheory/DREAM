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
    FVM::RadialGrid *rgrid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
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
 */
int_t EquationSystem::SetUnknown(const string& name, enum compregion region) {
    // TODO
    return -1;
}

