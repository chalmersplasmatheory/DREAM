/**
 * Process settings and build a Simulation object.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;

/**
 * Process the given settings and construct a
 * simulation object.
 *
 * s: Settings specifying how to construct the simulation.
 */
Simulation *SimulationGenerator::ProcessSettings(Settings *s) {
    // Construct grids
    enum OptionConstants::momentumgrid_type ht_type, re_type;
    FVM::Grid *fluidGrid   = ConstructRadialGrid(s);
    FVM::Grid *hottailGrid = ConstructHotTailGrid(s, fluidGrid->GetRadialGrid(), &ht_type);
    FVM::Grid *runawayGrid = ConstructRunawayGrid(s, fluidGrid->GetRadialGrid(), hottailGrid, &re_type);

    // Construct equation system
    EquationSystem *eqsys = ConstructEquationSystem(s, fluidGrid, ht_type, hottailGrid, re_type, runawayGrid);

    // Set up simulation
    Simulation *sim = new Simulation();
    sim->SetEquationSystem(eqsys);

    return sim;
}

/**
 * Convert an 'enum OptionConstants::momentumgrid_type' to
 * an 'enum Interpolator3D::momentumgrid_type'.
 */
enum FVM::Interpolator3D::momentumgrid_type SimulationGenerator::GetInterp3DMomentumGridType(
    enum OptionConstants::momentumgrid_type type
) {
    switch (type) {
        case OptionConstants::MOMENTUMGRID_TYPE_PXI: return FVM::Interpolator3D::GRID_PXI;
        case OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP: return FVM::Interpolator3D::GRID_PPARPPERP;

        default:
            throw SettingsException(
                "%s:%d: Asked to translate unrecognized momentum grid type: %d.",
                __FILE__, __LINE__, type
            );
    }
}
