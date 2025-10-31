/**
 * Implementation of class which is responsible for rebuilding collision quantities.
 */

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include <cmath>

using namespace DREAM;

/** 
 * Constructor
 */ 
CollisionQuantityHandler::CollisionQuantityHandler(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *u, 
    IonHandler *ih,  enum OptionConstants::momentumgrid_type gridtype,
	CollisionQuantity::collqty_settings *cqset
) : unknowns(u), ionHandler(ih), collQtySettings(cqset) {

	const len_t id_Tcold = u->GetUnknownID(OptionConstants::UQTY_T_COLD);
	const len_t id_ncold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);

    lnLambdaEE = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, collQtySettings, CollisionQuantity::LNLAMBDATYPE_EE, id_Tcold, id_ncold);
    lnLambdaEI = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, collQtySettings, CollisionQuantity::LNLAMBDATYPE_EI, id_Tcold, id_ncold);
    nuS   = new SlowingDownFrequency(grid, unknowns, ionHandler, lnLambdaEE,lnLambdaEI,gridtype, collQtySettings, id_Tcold, id_ncold);
    nuD   = new PitchScatterFrequency(grid, unknowns, ionHandler, lnLambdaEI,lnLambdaEE,gridtype, collQtySettings, id_Tcold, id_ncold);
    nuPar = new ParallelDiffusionFrequency(grid, unknowns, ionHandler, nuS,lnLambdaEE, gridtype, collQtySettings, id_Tcold, id_ncold);

	if (u->HasUnknown(OptionConstants::UQTY_T_HOT)) {
		const len_t id_Thot = u->GetUnknownID(OptionConstants::UQTY_T_HOT);
		const len_t id_nhot = u->GetUnknownID(OptionConstants::UQTY_N_HOT);

		lnLambdaEE = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, collQtySettings, CollisionQuantity::LNLAMBDATYPE_EE, id_Thot, id_nhot);
	}
}

/**
 * Destructor.
 */
CollisionQuantityHandler::~CollisionQuantityHandler(){
    delete nuS;
    delete nuD;
    delete nuPar;
    delete lnLambdaEE;
    delete lnLambdaEI;
	delete collQtySettings;
}

/**
 * Rebuilds all collision quantities
 */
void CollisionQuantityHandler::Rebuild() {
    lnLambdaEE->Rebuild();
    lnLambdaEI->Rebuild();
    nuS->Rebuild();
    nuD->Rebuild();
    nuPar->Rebuild();
}

/**
 * Notifies all collision quantities that the grid has been rebuilt.
 */
void CollisionQuantityHandler::gridRebuilt(){
    lnLambdaEE->GridRebuilt();
    lnLambdaEI->GridRebuilt();
    nuS->GridRebuilt();
    nuD->GridRebuilt();
    nuPar->GridRebuilt();
}
