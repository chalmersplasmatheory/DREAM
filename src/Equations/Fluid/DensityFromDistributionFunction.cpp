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
    const len_t N = fGrid->GetNCells();
    for (len_t i = 0; i < N; i++)
        this->integrand[i] = 1;
}


/**
 * Destructor.
 */
DensityFromDistributionFunction::~DensityFromDistributionFunction() { }

