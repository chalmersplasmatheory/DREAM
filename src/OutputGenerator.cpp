/**
 * General output generator interface.
 */

#include "DREAM/OutputGenerator.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
OutputGenerator::OutputGenerator(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns,
	IonHandler *ions, OtherQuantityHandler *oqty,
	EquationSystem *eqsys
) : grid(grid), unknowns(unknowns), ions(ions), oqty(oqty), eqsys(eqsys) { }

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

