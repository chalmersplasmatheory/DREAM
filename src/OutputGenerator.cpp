/**
 * General output generator interface.
 */

#include "DREAM/OutputGenerator.hpp"
#include "DREAM/Settings/Settings.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
OutputGenerator::OutputGenerator(
	EquationSystem *eqsys, bool savesettings
) : scalarGrid(eqsys->GetScalarGrid()), fluidGrid(eqsys->GetFluidGrid()),
    hottailGrid(eqsys->GetHotTailGrid()), runawayGrid(eqsys->GetRunawayGrid()),
    unknowns(eqsys->GetUnknownHandler()), ions(eqsys->GetIonHandler()),
    oqty(eqsys->GetOtherQuantityHandler()), eqsys(eqsys),
    savesettings(savesettings) { }

OutputGenerator::~OutputGenerator() {}

/**
 * Generate the desired output.
 *
 * current: If true, saves only unknown data from the current
 *          iteration or time step.
 */
void OutputGenerator::Save(bool current) {
    // Save grids
    this->SaveGrids("grid", current);

    // Save unknowns
    this->SaveUnknowns("eqsys", current);

    // Save ion metadata
    this->SaveIonMetaData("ionmeta");

    // Save settings
    if (this->savesettings)
        this->SaveSettings("settings");

    // Save solver statistics
    this->SaveSolverData("solver");

    // Save timing information
    this->SaveTimings("timings");

    // Save "other" quantities (if any)
    if (oqty->GetNRegistered() > 0 && !current)
        this->SaveOtherQuantities("other");
	
	// Save information about the code
	this->SaveCodeInfo("code");
}

