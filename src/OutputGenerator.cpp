/**
 * General output generator interface.
 */

#include "DREAM/OutputGenerator.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
OutputGenerator::OutputGenerator(
	EquationSystem *eqsys
) : scalarGrid(eqsys->GetScalarGrid()), fluidGrid(eqsys->GetFluidGrid()),
    hottailGrid(eqsys->GetHotTailGrid()), runawayGrid(eqsys->GetRunawayGrid()),
    unknowns(eqsys->GetUnknownHandler()), ions(eqsys->GetIonHandler()),
    oqty(eqsys->GetOtherQuantityHandler()), eqsys(eqsys) { }

OutputGenerator::~OutputGenerator() {}

/**
 * Generate the desired output.
 *
 * current: If true, saves only unknown data from the current
 *          iteration or time step.
 */
void OutputGenerator::Save(bool current) {
    // TODO Save settings

    // Save grids
    this->SaveGrids("grid", current);
    
    // Save unknowns
    this->SaveUnknowns("eqsys", current);

    // Save ion metadata
    this->SaveIonMetaData("ionmeta");

    // Save timing information
    this->SaveTimings("timings");

    // Save "other" quantities (if any)
    if (oqty->GetNRegistered() > 0 && !current)
        this->SaveOtherQuantities("other");
}

