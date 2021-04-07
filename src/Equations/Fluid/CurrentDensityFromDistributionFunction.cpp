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

    SetName("CurrentDensityFromDistributionFunction");

    this->scaleFactor = scaleFactor;
    // Build moment integrand
    this->GridRebuilt();
}

/**
 * Method that is called whenever the grid is rebuilt. 
 */
bool CurrentDensityFromDistributionFunction::GridRebuilt() {
    if (this->MomentQuantity::GridRebuilt()) {
        // XXX assumes same momentumgrid at all radii.
        FVM::MomentumGrid *mg;
        FVM::RadialGrid *rGrid = fGrid->GetRadialGrid();

        len_t np1, np2, ind;
        len_t offset = 0;
        for(len_t ir = 0; ir<rGrid->GetNr(); ir++){
            mg = fGrid->GetMomentumGrid(ir);
            np1 = mg->GetNp1();
            np2 = mg->GetNp2();
            const real_t *Vp = fGrid->GetVp(ir);
            const real_t VpVol = fGrid->GetVpVol(ir);
            real_t xi0Trapped = rGrid->GetXi0TrappedBoundary(ir);
            for(len_t ip1 = 0; ip1<np1; ip1++){
                for(len_t ip2 = 0; ip2<np2; ip2++){
                    ind = ip2*np1+ip1;
                    real_t p = mg->GetP(ip1,ip2);
                    real_t v = Constants::c * p / mg->GetGamma(ip1,ip2);
                    
                    real_t xi1 = mg->GetXi0_f2(ip1,ip2);
                    real_t xi2 = mg->GetXi0_f2(ip1,ip2+1);
                    if(xi1>xi2){ // if xi's are decreasing, switch upper and lower
                        real_t xi_tmp = xi1;
                        xi1 = xi2;
                        xi2 = xi_tmp;
                    }
                    // xi0Factor is the cell average of xi0 over the passing region
                    real_t xi0Factor = 0;
                    // cell entirely in passing region or entire trapped region in cell
                    if( xi2<=-xi0Trapped || xi1>=xi0Trapped || (xi1<=-xi0Trapped && xi2>=xi0Trapped))  
                        xi0Factor = mg->GetXi0(ip1,ip2);

                    // replace Vp/VpVol in MomentQuantity by 2pi*p^2 jacobian
                    real_t Jacobian = 0;
                    if(xi0Factor) 
                        Jacobian = 2*M_PI * p*p * VpVol / Vp[ind]; 
                    this->integrand[offset+ind] = Jacobian * Constants::ec * v * xi0Factor;
                }
            }
            offset += np1*np2;
        }

        return true;
    } else return false;
}
