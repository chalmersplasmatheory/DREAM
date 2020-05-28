/**
 * Implementation of collision-rate calculator that calculates
 * various collision rates and related quantities, such as runaway growth rates.
*/



#include "DREAM/Equations/CollisionQuantityHandler.hpp"
//#include "DREAM/Constants.hpp"
//#include "DREAM/Settings/OptionConstants.hpp"
//#include "FVM/UnknownQuantityHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include <cmath>
//#include <gsl/gsl_roots.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_min.h>
//#include <iostream>
//#include <fstream>


using namespace DREAM;

const len_t  CollisionQuantityHandler::conductivityLenT = 14;
const len_t  CollisionQuantityHandler::conductivityLenZ = 6;
const real_t CollisionQuantityHandler::conductivityBraams[conductivityLenZ*conductivityLenT] = {3.75994, 3.7549, 3.7492, 3.72852, 3.6842, 3.57129, 3.18206, 2.65006, 2.03127, 1.33009, 0.94648, 0.67042, 0.42422, 0.29999, 7.42898, 7.27359, 7.12772, 6.73805, 6.20946, 5.43667, 4.13733, 3.13472, 2.27862, 1.45375, 1.02875, 0.72743, 0.46003, 0.32528, 8.7546, 8.53281, 8.32655, 7.78445, 7.06892, 6.06243, 4.47244, 3.32611, 2.39205, 1.51805, 1.07308, 0.75853, 0.47965, 0.33915, 10.39122, 10.07781, 9.78962, 9.04621, 8.09361, 6.80431, 4.8805, 3.57303, 2.54842, 1.61157, 1.13856, 0.80472, 0.50885, 0.35979, 11.33006, 10.95869, 10.61952, 9.75405, 8.66306, 7.21564, 5.11377, 3.72206, 2.64827, 1.67382, 1.18263, 0.83593, 0.52861, 0.37377, 12.76615, 12.29716, 11.87371, 10.81201, 9.50746, 7.82693, 5.47602, 3.96944, 2.82473, 1.7887, 1.2649, 0.89443, 0.56569, 0.4};
const real_t CollisionQuantityHandler::conductivityTmc2[conductivityLenT] = {0,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100};
const real_t CollisionQuantityHandler::conductivityX[conductivityLenZ]    = {0,0.090909090909091,0.166666666666667,0.333333333333333,0.5,1};

/** 
 * Constructor
 */ 
CollisionQuantityHandler::CollisionQuantityHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantity::collqty_settings *cqset){
    ionHandler = ih;
    grid       = g;
    unknowns   = u;
    settings   = cqset;
    gridtype   = mgtype;

    this->nr   = g->GetNr();

    lnLambdaEE = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, settings,CollisionQuantity::LNLAMBDATYPE_EE);
    lnLambdaEI = new CoulombLogarithm(grid, unknowns, ionHandler, gridtype, settings,CollisionQuantity::LNLAMBDATYPE_EI);
    nuS   = new SlowingDownFrequency(grid, unknowns, ionHandler, lnLambdaEE,lnLambdaEI,gridtype, settings);
    nuD   = new PitchScatterFrequency(grid, unknowns, ionHandler, lnLambdaEI,lnLambdaEE,gridtype, settings);
    nuPar = new ParallelDiffusionFrequency(grid, unknowns, ionHandler, nuS,lnLambdaEE, gridtype, settings);
    REFluid = new RunawayFluid(grid,unknowns,nuS,nuD,lnLambdaEE,settings);


    const gsl_interp2d_type *gsl_T = gsl_interp2d_bilinear; 
    gsl_cond = gsl_interp2d_alloc(gsl_T, conductivityLenT,conductivityLenZ);
    gsl_xacc = gsl_interp_accel_alloc();
    gsl_yacc = gsl_interp_accel_alloc();

    gsl_interp2d_init(gsl_cond, conductivityTmc2, conductivityX, conductivityBraams,conductivityLenT,conductivityLenZ);
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
    delete [] REFluid;

    gsl_interp2d_free(gsl_cond);
    gsl_interp_accel_free(gsl_xacc);
    gsl_interp_accel_free(gsl_yacc);
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
    bool useApproximateEceffCalculation = false;
    REFluid->Rebuild(useApproximateEceffCalculation);

    // for testing:
    // RunawayFluid *REFluid2 = new RunawayFluid(grid,unknowns,nuS,nuD,lnLambdaEE,settings);
    // REFluid2->Rebuild(true);

}

void CollisionQuantityHandler::gridRebuilt(){
    lnLambdaEE->GridRebuilt();
    lnLambdaEI->GridRebuilt();
    nuS->GridRebuilt();
    nuD->GridRebuilt();
    nuPar->GridRebuilt();
}

real_t CollisionQuantityHandler::evaluateElectricalConductivity(len_t ir){
    len_t id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    const real_t T_SI = T_cold[ir] * Constants::ec;
    const real_t *Zeff = ionHandler->evaluateZeff();

    real_t sigmaBar = gsl_interp2d_eval(gsl_cond, conductivityTmc2, conductivityX, conductivityBraams, 
                T_SI / (Constants::me * Constants::c * Constants::c), 1/(1+Zeff[ir]), gsl_xacc, gsl_yacc  );
    
    real_t BraamsConductivity = 4*M_PI*Constants::eps0*Constants::eps0 * T_SI*sqrt(T_SI) / 
            (sqrt(Constants::me) * Constants::ec * Constants::ec * lnLambdaEE->GetLnLambdaT(ir) ) * sigmaBar;
    delete [] Zeff;
    return BraamsConductivity;
}

