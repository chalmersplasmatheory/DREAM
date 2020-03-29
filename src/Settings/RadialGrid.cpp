/**
 * Construction of the radial grid.
 */

#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"


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
    s->DefineSetting(RADIALGRID "/type", "Type of radial grid", (int_t)RADIALGRID_TYPE_CYLINDRICAL);

    // CylindricalRadialGrid
    s->DefineSetting(RADIALGRID "/a",  "Tokamak minor radius", (real_t)0.5);
    s->DefineSetting(RADIALGRID "/B0", "On-axis magnetic field strength", (real_t)1.0);
    s->DefineSetting(RADIALGRID "/r0", "Inner-most radius to simulate (on flux-grid)", (real_t)0.0);
}

/**
 * Construct a radial grid according to the
 * given specification.
 *
 * s: Settings object specifying how to construct
 *    the radial grid.
 */
FVM::RadialGrid *SimulationGenerator::ConstructRadialGrid(Settings *s) {
    enum radialgrid_type type = (enum radialgrid_type)s->GetInteger(RADIALGRID "/type");
    int_t nr = s->GetInteger(RADIALGRID "/nr");

    switch (type) {
        case RADIALGRID_TYPE_CYLINDRICAL:
            return ConstructRadialGrid_Cylindrical(nr, s);

        default:
            throw SettingsException(
                "Unrecognized radial grid type specified: " INT_T_PRINTF_FMT ".",
                type
            );
    }
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

