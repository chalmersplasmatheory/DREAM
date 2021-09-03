/**
 * Implementation of an operator which evaluates the total heat flux from a distribution
 */

#include "DREAM/Equations/Fluid/HeatFluxFromDistributionFunction.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
HeatFluxFromDistributionFunction::HeatFluxFromDistributionFunction(
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid, len_t id_n, len_t id_f,
    FVM::UnknownQuantityHandler *u, real_t pThreshold, pThresholdMode pMode, real_t scaleFactor)
     : MomentQuantity(densityGrid, distributionGrid, id_n, id_f, u, pThreshold, pMode) {

    this->scaleFactor = scaleFactor;
    // Build moment integrand
    this->GridRebuilt();
}

/**
 * Method that is called whenever the grid is rebuilt. 
 */
bool HeatFluxFromDistributionFunction::GridRebuilt() {
    if (this->MomentQuantity::GridRebuilt()) {
        FVM::MomentumGrid *mg;
        FVM::RadialGrid *rGrid = fGrid->GetRadialGrid();

        len_t np1, np2, ind;
        len_t offset = 0;
        for(len_t ir = 0; ir<rGrid->GetNr(); ir++){
            mg = fGrid->GetMomentumGrid(ir);
            np1 = mg->GetNp1();
            np2 = mg->GetNp2();
            for(len_t ip1 = 0; ip1<np1; ip1++){
                for(len_t ip2 = 0; ip2<np2; ip2++){
                    ind = ip2*np1+ip1;
                    real_t p = mg->GetP(ip1,ip2);
                    real_t v = Constants::c * p / mg->GetGamma(ip1,ip2);
                    this->integrand[offset+ind] = scaleFactor * Constants::me*Constants::c*Constants::c*(mg->GetGamma(ip1,ip2)-1) * v;
                }
            }
            offset += np1*np2;
        }

        return true;
    } else return false;
}
