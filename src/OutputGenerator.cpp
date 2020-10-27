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
 */
void OutputGenerator::Save() {
    // TODO Save settings

    // Save grids
    this->SaveGrids("grid");
    
    // Save unknowns
    this->SaveUnknowns("eqsys");

    // Save ion metadata
    this->SaveIonMetaData("ionmeta");

    // Save timing information
    this->SaveTimings("timings");

    // Save "other" quantities (if any)
    if (oqty->GetNRegistered() > 0)
        this->SaveOtherQuantities("other");
}

