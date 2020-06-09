/**
 * Implementation of class which is responsible for rebuilding collision quantities.
 */

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include <cmath>

using namespace DREAM;

/** 
 * Constructor
 */ 
CollisionQuantityHandler::CollisionQuantityHandler(FVM::Grid *grid, FVM::UnknownQuantityHandler *u, IonHandler *ih,  enum OptionConstants::momentumgrid_type gridtype,  struct CollisionQuantity::collqty_settings *cqset){
    ionHandler = ih;
    unknowns   = u;

    lnLambdaEE = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, cqset,CollisionQuantity::LNLAMBDATYPE_EE);
    lnLambdaEI = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, cqset,CollisionQuantity::LNLAMBDATYPE_EI);
    nuS   = new SlowingDownFrequency(grid, unknowns, ionHandler, lnLambdaEE,lnLambdaEI,gridtype, cqset);
    nuD   = new PitchScatterFrequency(grid, unknowns, ionHandler, lnLambdaEI,lnLambdaEE,gridtype, cqset);
    nuPar = new ParallelDiffusionFrequency(grid, unknowns, ionHandler, nuS,lnLambdaEE, gridtype, cqset);

    }

/**
 * Destructor.
 */
CollisionQuantityHandler::~CollisionQuantityHandler(){
    delete [] nuS;
    delete [] nuD;
    delete [] nuPar;
    delete [] lnLambdaEE;
    delete [] lnLambdaEI;

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
