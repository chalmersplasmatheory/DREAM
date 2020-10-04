/**
 * Construction of the radial grid.
 */

#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "FVM/Grid/EmptyRadialGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace DREAM;

// Module name (to give compile-time error if misspelled,
// instead of a run-time error)
#define RADIALGRID "radialgrid"

/**
 * Define options for the radial simulation grid.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_RadialGrid(Settings *s) {

    // Radial grid
    s->DefineSetting(RADIALGRID "/nr",   "Number of radial (distribution) grid points", (int_t)1, true);
    s->DefineSetting(RADIALGRID "/type", "Type of radial grid", (int_t)OptionConstants::RADIALGRID_TYPE_CYLINDRICAL);

    // CylindricalRadialGrid
    s->DefineSetting(RADIALGRID "/a",  "Tokamak minor radius", (real_t)0.5);
    s->DefineSetting(RADIALGRID "/B0", "On-axis magnetic field strength", (real_t)1.0);
    s->DefineSetting(RADIALGRID "/r0", "Inner-most radius to simulate (on flux-grid)", (real_t)0.0);

    // AnalyticBRadialGridGenerator
    s->DefineSetting(RADIALGRID "/R0", "Tokamak major radius", (real_t)2.0);

    DefineDataR(RADIALGRID, s, "delta");
    DefineDataR(RADIALGRID, s, "kappa");

    // TODO other equilibrium parameters?
}

/**
 * Construct a radial grid according to the
 * given specification.
 *
 * s: Settings object specifying how to construct
 *    the radial grid.
 */
FVM::Grid *SimulationGenerator::ConstructRadialGrid(Settings *s) {
    enum OptionConstants::radialgrid_type type = (enum OptionConstants::radialgrid_type)s->GetInteger(RADIALGRID "/type");
    int_t nr = s->GetInteger(RADIALGRID "/nr");

    FVM::RadialGrid *rg;
    switch (type) {
        case OptionConstants::RADIALGRID_TYPE_CYLINDRICAL:
            rg = ConstructRadialGrid_Cylindrical(nr, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized radial grid type specified: " INT_T_PRINTF_FMT ".",
                type
            );
    }
    
    return new FVM::Grid(rg, new FVM::EmptyMomentumGrid(rg));
}



/**
 * Construct a radial grid according to the
 * given specification.
 *
 * s: Settings object specifying how to construct
 *    the radial grid.
 */
FVM::Grid *SimulationGenerator::ConstructScalarGrid() {
//    auto *emptyGridGenerator = new FVM::EmptyRadialGridGenerator();
//   FVM::RadialGrid *rg = new FVM::RadialGrid(emptyGridGenerator);
    FVM::RadialGrid *rg = new FVM::EmptyRadialGrid();
    return new FVM::Grid(rg, new FVM::EmptyMomentumGrid(rg));
}



/*************************************
 * SPECIFIC RADIAL GRID CONSTRUCTION *
 *************************************/
/**
 * Construct a radial grid that approximates the tokamak
 * as a straight cylinder.
 * 
 * nr: Number of radial (distribution) grid points.
 * s:  Settings object specifying how to construct the grid.
 */
FVM::RadialGrid *SimulationGenerator::ConstructRadialGrid_Cylindrical(const int_t nr, Settings *s) {
    real_t a  = s->GetReal(RADIALGRID "/a");
    real_t B0 = s->GetReal(RADIALGRID "/B0");
    real_t r0 = s->GetReal(RADIALGRID "/r0");

    auto *crgg = new FVM::CylindricalRadialGridGenerator(nr, B0, r0, a);
    return new FVM::RadialGrid(crgg);
}

