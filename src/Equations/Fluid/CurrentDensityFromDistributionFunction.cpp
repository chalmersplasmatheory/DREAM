/**
 * Implementation of an operator which evaluates the parallel current density
 * moment j_||/(B/Bmin) normalized to the local magnetic field strength
 * of a given distribution function.
 */

#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
CurrentDensityFromDistributionFunction::CurrentDensityFromDistributionFunction(
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid, len_t id_n, len_t id_f
) : MomentQuantity(densityGrid, distributionGrid, id_n, id_f) {
    
    // Build moment integrand
    this->GridRebuilt();
}


/**
 * Destructor.
 */
CurrentDensityFromDistributionFunction::~CurrentDensityFromDistributionFunction() { }


/**
 * Method that is called whenever the grid is rebuilt. We only
 * need to rebuild this EquationTerm if the total number of grid
 * cells changes.
 */
bool CurrentDensityFromDistributionFunction::GridRebuilt() {
    if (this->MomentQuantity::GridRebuilt()) {
        // XXX assumes same momentumgrid at all radii.
        FVM::MomentumGrid *mg;
        FVM::RadialGrid *rGrid = fGrid->GetRadialGrid();

        len_t np1, np2, ind;
        real_t v, xi0, geometricFactor;
        len_t offset = 0;
        for(len_t ir = 0; ir<rGrid->GetNr(); ir++){
            mg = fGrid->GetMomentumGrid(ir);
            np1 = mg->GetNp1();
            np2 = mg->GetNp2();
            for(len_t ip1 = 0; ip1<np1; ip1++){
                for(len_t ip2 = 0; ip2<np2; ip2++){
                    ind = offset+ip2*np1+ip1;
                    v = Constants::c *mg->GetP(ip1,ip2)/mg->GetGamma(ip1,ip2);
                    xi0 = mg->GetXi0(ip1,ip2);
                    geometricFactor = grid->IsTrapped(ir,ip1,ip2);
                    this->integrand[ind] = Constants::ec * v * xi0 * geometricFactor;
                }
            }
            offset += np1*np2;
        }

        return true;
    } else return false;
}
