/**
 * Process settings and build a Simulation object.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace DREAM;

/**
 * Process the given settings and construct a
 * simulation object.
 *
 * s: Settings specifying how to construct the simulation.
 */
Simulation *SimulationGenerator::ProcessSettings(Settings *s) {
    // Construct grids
    FVM::Grid *fluidGrid   = ConstructRadialGrid(s);
    FVM::Grid *hottailGrid = ConstructHotTailGrid(s, fluidGrid->GetRadialGrid());
    FVM::Grid *runawayGrid = ConstructRunawayGrid(s, fluidGrid->GetRadialGrid(), hottailGrid);

    // Construct equation system
    EquationSystem *eqsys = ConstructEquationSystem(s, fluidGrid, hottailGrid, runawayGrid);

    // Set up simulation
    Simulation *sim = new Simulation();
    sim->SetEquationSystem(eqsys);

    return sim;
}

