/**
 * Implementation of an operator which evaluates the density
 * moment of a given distribution function.
 */

#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
DensityFromDistributionFunction::DensityFromDistributionFunction(
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid, len_t id_n, len_t id_f,
    FVM::UnknownQuantityHandler *u, real_t pThreshold, pThresholdMode pMode,
    xiIntegralMode xiMode
) : MomentQuantity(densityGrid, distributionGrid, id_n, id_f, u, pThreshold, pMode, xiMode) {

    SetName("DensityFromDistributionFunction");

    // Build moment integrand
    this->GridRebuilt();
}

/**
 * Method that is called whenever the grid is rebuilt. We only
 * need to rebuild this EquationTerm if the total number of grid
 * cells changes.
 */
bool DensityFromDistributionFunction::GridRebuilt() {
    if (this->MomentQuantity::GridRebuilt()) {
        for (len_t i = 0; i < this->nIntegrand; i++)
            this->integrand[i] = 1;
        return true;
    } else 
        return false;
}
