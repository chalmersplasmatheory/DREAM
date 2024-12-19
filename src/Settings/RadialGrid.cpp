/**
 * Construction of the radial grid.
 */

#include <string>
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/NumericBRadialGridGenerator.hpp"
#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "FVM/Grid/EmptyRadialGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace DREAM;
using namespace std;

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

    s->DefineSetting(RADIALGRID "/a",  "Tokamak minor radius", (real_t)0.5);
    s->DefineSetting(RADIALGRID "/r0", "Inner-most radius to simulate (on flux-grid)", (real_t)0.0);
    s->DefineSetting(RADIALGRID "/wall_radius",  "Tokamak wall minor radius", (real_t)0.5);

	s->DefineSetting(RADIALGRID "/custom_grid", "Flag indicating whether to use a custom radial grid or not", (bool)false);
    s->DefineSetting(RADIALGRID "/r_f", "Grid points of the radial flux grid", 0, (real_t*) nullptr);

    // CylindricalRadialGrid
    s->DefineSetting(RADIALGRID "/B0", "On-axis magnetic field strength", (real_t)1.0);

    // AnalyticBRadialGridGenerator
    s->DefineSetting(RADIALGRID "/R0", "Tokamak major radius", (real_t)2.0);
    s->DefineSetting(RADIALGRID "/ntheta", "Number of poloidal angle grid points to use for bounce averages", (int_t)30);
	s->DefineSetting(RADIALGRID "/ntheta_out", "Number of poloidal angle grid points to use on output flux surface grid", (int_t)120);

    DefineDataR(RADIALGRID, s, "delta");    // Triangularity
    DefineDataR(RADIALGRID, s, "Delta");    // Shafranov shift
    DefineDataR(RADIALGRID, s, "kappa");    // Elongation
    DefineDataR(RADIALGRID, s, "GOverR0");        // G/R0 = (R/R0)*Bphi
    DefineDataR(RADIALGRID, s, "psi_p0");   // Reference poloidal flux (normalized to R0)
    // Magnetic ripple effects
    DefineOptions_f_ripple(RADIALGRID, s);
	// Time-varying B operator
	DefineOptions_f_timevaryingb(RADIALGRID, s);

    // NumericBRadialGridGenerator
    s->DefineSetting(RADIALGRID "/filename", "Name of file containing the magnetic field data", (string)"");
    s->DefineSetting(RADIALGRID "/fileformat", "Format used for storing the magnetic field data", (int_t)OptionConstants::RADIALGRID_NUMERIC_FORMAT_LUKE);
}

/**
 * Define options for the magnetic ripple modelling.
 */
void SimulationGenerator::DefineOptions_f_ripple(const string& mod, Settings *s) {
    s->DefineSetting(mod + "/ripple/ncoils", "Number of toroidal magnetic field coils", (int_t)0);
    s->DefineSetting(mod + "/ripple/deltacoils", "Distance between magnetic field coils (alternative to ncoils)", (real_t)0);
    
    s->DefineSetting(mod + "/ripple/m", "Poloidal mode numbers", 0, (int_t*)nullptr);
    s->DefineSetting(mod + "/ripple/n", "Toroidal mode numbers", 0, (int_t*)nullptr);

    // Define perturbation data
    DefineDataIonRT(mod, s, "ripple");
}

/**
 * Define options for the time-varying magnetic field operator.
 */
void SimulationGenerator::DefineOptions_f_timevaryingb(
	const string& mod, Settings *s
) {
	DefineDataT(mod, s, "dlnB0dt");
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

        case OptionConstants::RADIALGRID_TYPE_TOROIDAL_ANALYTICAL:
            rg = ConstructRadialGrid_ToroidalAnalytical(nr, s);
            break;

        case OptionConstants::RADIALGRID_TYPE_NUMERICAL:
            rg = ConstructRadialGrid_Numerical(nr, s);
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
    real_t B0 = s->GetReal(RADIALGRID "/B0");
	len_t ntheta_out = s->GetInteger(RADIALGRID "/ntheta_out");
	bool custom_grid = s->GetBool(RADIALGRID "/custom_grid");

    FVM::CylindricalRadialGridGenerator *crgg;
    if(!custom_grid){
        real_t a  = s->GetReal(RADIALGRID "/a");
        real_t r0 = s->GetReal(RADIALGRID "/r0");
        crgg = new FVM::CylindricalRadialGridGenerator(nr, B0, r0, a, ntheta_out);
    } else {
        len_t len_rf; // equals nr+1 of the simulation
        const real_t *r_f = s->GetRealArray(RADIALGRID "/r_f", 1, &len_rf);
        crgg = new FVM::CylindricalRadialGridGenerator(r_f, len_rf-1, B0, ntheta_out);
    }
    return new FVM::RadialGrid(crgg);
}

/**
 * Construct a toroidal radial grid using an analytical
 * model for the magnetic field.
 *
 * nr: Number of radial (distribution) grid points.
 * s:  Settings object specifying how to construct the grid.
 */
FVM::RadialGrid *SimulationGenerator::ConstructRadialGrid_ToroidalAnalytical(const int_t nr, Settings *s) {
    real_t R0 = s->GetReal(RADIALGRID "/R0");
    len_t ntheta_interp = s->GetInteger(RADIALGRID "/ntheta");
	len_t ntheta_out = s->GetInteger(RADIALGRID "/ntheta_out");

    FVM::AnalyticBRadialGridGenerator::shape_profiles *shapes =
        new FVM::AnalyticBRadialGridGenerator::shape_profiles;

    shapes->GOverR0 = s->GetRealArray(RADIALGRID "/GOverR0/x", 1, &shapes->nG);
    shapes->G_r     = s->GetRealArray(RADIALGRID "/GOverR0/r", 1, &shapes->nG);
    shapes->delta   = s->GetRealArray(RADIALGRID "/delta/x", 1, &shapes->ndelta);
    shapes->delta_r = s->GetRealArray(RADIALGRID "/delta/r", 1, &shapes->ndelta);
    shapes->Delta   = s->GetRealArray(RADIALGRID "/Delta/x", 1, &shapes->nDelta);
    shapes->Delta_r = s->GetRealArray(RADIALGRID "/Delta/r", 1, &shapes->nDelta);
    shapes->kappa   = s->GetRealArray(RADIALGRID "/kappa/x", 1, &shapes->nkappa);
    shapes->kappa_r = s->GetRealArray(RADIALGRID "/kappa/r", 1, &shapes->nkappa);
    shapes->psi     = s->GetRealArray(RADIALGRID "/psi_p0/x", 1, &shapes->npsi);
    shapes->psi_r   = s->GetRealArray(RADIALGRID "/psi_p0/r", 1, &shapes->npsi);

    FVM::AnalyticBRadialGridGenerator*abrg;
    if(nr!=0){
        real_t a  = s->GetReal(RADIALGRID "/a");
        real_t r0 = s->GetReal(RADIALGRID "/r0");
        abrg = new FVM::AnalyticBRadialGridGenerator(
            nr, r0, a, R0, ntheta_interp, shapes, ntheta_out
        );
    } else {
        len_t len_rf; // equals nr+1 of the simulation
        const real_t *r_f = s->GetRealArray(RADIALGRID "/r_f", 1, &len_rf);
        abrg = new FVM::AnalyticBRadialGridGenerator(
            r_f, len_rf-1, R0, ntheta_interp, shapes, ntheta_out
        );
    }

    return new FVM::RadialGrid(abrg);
}

/**
 * Construct a radial grid based on a numeric magnetic field.
 *
 * nr: Number of radial (distribution) grid points.
 * s:  Settings object specifying how to construct the grid.
 */
FVM::RadialGrid *SimulationGenerator::ConstructRadialGrid_Numerical(
    const int_t nr, Settings *s
) {
    len_t ntheta_interp = s->GetInteger(RADIALGRID "/ntheta");

    const string filename = s->GetString(RADIALGRID "/filename");
    enum OptionConstants::radialgrid_numeric_format frmt =
        (enum OptionConstants::radialgrid_numeric_format)s->GetInteger(RADIALGRID "/fileformat");

    enum FVM::NumericBRadialGridGenerator::file_format nbrg_frmt;
    switch (frmt) {
        case OptionConstants::RADIALGRID_NUMERIC_FORMAT_LUKE:
            nbrg_frmt = FVM::NumericBRadialGridGenerator::FILE_FORMAT_LUKE;
            break;

        default:
            throw SettingsException(
                "NumericBRadialGrid settings: Unrecognized magnetic field data format: %d.",
                frmt
            );
    }

    FVM::NumericBRadialGridGenerator *nbrg;

    // Uniform radial grid
    if (nr != 0) {
        real_t a  = s->GetReal(RADIALGRID "/a");
        real_t r0 = s->GetReal(RADIALGRID "/r0");

        nbrg = new FVM::NumericBRadialGridGenerator(
            nr, r0, a, filename, nbrg_frmt, ntheta_interp
        );

    // Custom radial grid
    } else {
        len_t len_rf; // equals nr+1 of the simulation
        const real_t *r_f = s->GetRealArray(RADIALGRID "/r_f", 1, &len_rf);

        nbrg = new FVM::NumericBRadialGridGenerator(
            r_f, len_rf-1, filename, nbrg_frmt, ntheta_interp
        );
    }

    return new FVM::RadialGrid(nbrg);
}

