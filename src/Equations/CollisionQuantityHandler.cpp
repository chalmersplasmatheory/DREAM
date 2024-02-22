/**
 * Implementation of class which is responsible for rebuilding collision quantities.
 */

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include <cmath>

using namespace DREAM;

/** 
 * Constructor
 */ 
CollisionQuantityHandler::CollisionQuantityHandler(FVM::Grid *grid, FVM::UnknownQuantityHandler *u, 
    IonHandler *ih,  enum OptionConstants::momentumgrid_type gridtype,  CollisionQuantity::collqty_settings *cqset)
        : unknowns(u), ionHandler(ih),collQtySettings(cqset) 
{
    lnLambdaEE = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, collQtySettings,CollisionQuantity::LNLAMBDATYPE_EE);
    lnLambdaEI = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, collQtySettings,CollisionQuantity::LNLAMBDATYPE_EI);
    nuS   = new SlowingDownFrequency(grid, unknowns, ionHandler, lnLambdaEE,lnLambdaEI,gridtype, collQtySettings);
    nuD   = new PitchScatterFrequency(grid, unknowns, ionHandler, lnLambdaEI,lnLambdaEE,gridtype, collQtySettings);
    nuPar = new ParallelDiffusionFrequency(grid, unknowns, ionHandler, nuS,lnLambdaEE, gridtype, collQtySettings);
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
