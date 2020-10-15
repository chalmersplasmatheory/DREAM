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
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid, len_t id_n, len_t id_f,
    FVM::UnknownQuantityHandler *u, real_t pThreshold, pThresholdMode pMode, real_t scaleFactor)
     : MomentQuantity(densityGrid, distributionGrid, id_n, id_f, u, pThreshold, pMode) {

    this->scaleFactor = scaleFactor;
    // Build moment integrand
    this->GridRebuilt();
}

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
            real_t xi0Trapped = rGrid->GetXi0TrappedBoundary(ir);
            for(len_t ip1 = 0; ip1<np1; ip1++){
                for(len_t ip2 = 0; ip2<np2; ip2++){
                    ind = offset+ip2*np1+ip1;
                    v = Constants::c *mg->GetP(ip1,ip2)/mg->GetGamma(ip1,ip2);
                    xi0 = mg->GetXi0(ip1,ip2);

                    // the geometricFactor is the fraction of the cell that lies in the passing region
                    // This is a compacted method of evaluating the cell-averaged (over pitch) bounce integral
                    // XXX: it assumes p-xi grid to work optimally (where _f2 is the pitch flux grid)
                    geometricFactor = 1;
                    if(xi0Trapped){ 
                        real_t xi1 = mg->GetXi0_f2(ip1,ip2);
                        real_t xi2 = mg->GetXi0_f2(ip1,ip2+1);
                        real_t dxiBarTrapped = std::min(xi2,xi0Trapped) - std::max(xi1,-xi0Trapped); // pitch interval that overlaps with trapped region
                        if(dxiBarTrapped>0)
                            geometricFactor = 1 - dxiBarTrapped / (xi2-xi1);
                    }
                    this->integrand[ind] = Constants::ec * v * xi0 * geometricFactor;
                }
            }
            offset += np1*np2;
        }

        return true;
    } else return false;
}
