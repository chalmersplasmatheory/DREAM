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
        real_t *const*bounceAverage = fGrid->GetBA_xiOverBR2();
        const real_t *fluxSurfaceAverage = rGrid->GetFSA_1OverR2();
        //for (len_t i = 0; i < this->nIntegrand; i++){ // i = j0*np1 + i0. j0 = i/np1. i0 = i%np1 
        for(len_t ir = 0; ir<rGrid->GetNr(); ir++){
            mg = fGrid->GetMomentumGrid(ir);
            np1 = mg->GetNp1();
            np2 = mg->GetNp2();
            for(len_t ip1 = 0; ip1<np1; ip1++)
                for(len_t ip2 = 0; ip2<np2; ip2++){
                    // the geometric factor equals 1 for passing particles and 0 for trapped particles. 
                    // It should be identical to rGrid->GetIsTrapped(...).
                    //if(IsTrapped(ir,ip1,ip2))
                    //    this->integrand[ind] = 0;
                    //else ...
                    ind = ir*np1*np2+ip2*np1+ip1;
                    v = Constants::c *mg->GetP(ip1,ip2)/mg->GetGamma(ip1,ip2);
                    xi0 = mg->GetXi0(ip1,ip2);
                    geometricFactor = bounceAverage[ir][ip2*np1+ip1] / fluxSurfaceAverage[ir]; 
                    this->integrand[ind] = scaleFactor * Constants::ec * v * xi0 * geometricFactor;
                }
        }

        return true;
    } else return false;
}
