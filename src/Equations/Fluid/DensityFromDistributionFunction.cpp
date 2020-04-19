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
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid, len_t id_n, len_t id_f
) : MomentQuantity(densityGrid, distributionGrid, id_n, id_f) {
    
    // Build moment integrand
    this->GridRebuilt();
}


/**
 * Destructor.
 */
DensityFromDistributionFunction::~DensityFromDistributionFunction() { }


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
    } else return false;
}
