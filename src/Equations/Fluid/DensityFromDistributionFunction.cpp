/**
 * Implementation of an operator which evaluates the density
 * moment of a given distribution function.
 */

#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
DensityFromDistributionFunction::DensityFromDistributionFunction(
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid, len_t id_n, len_t id_f,
    real_t pThreshold, pThresholdMode pMode
) : MomentQuantity(densityGrid, distributionGrid, id_n, id_f), pThreshold(pThreshold), pMode(pMode) 
{
    // Build moment integrand
    this->GridRebuilt();
}


/**
 * Destructor.
 */
DensityFromDistributionFunction::~DensityFromDistributionFunction() { }

/**
 * Rebuilds the integrand of this term.
 */
void DensityFromDistributionFunction::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler* unknowns) {
    ResetIntegrand();
    FVM::MomentumGrid *mg;
    const len_t nr = fGrid->GetNr();
    len_t np1, np2, ind;
    for(len_t ir = 0; ir<nr; ir++){
        mg = fGrid->GetMomentumGrid(ir);
        np1 = mg->GetNp1();
        np2 = mg->GetNp2();
        real_t p0 = GetThreshold(ir,unknowns);
        len_t iMin = (p0>0); // XXX: assumes p-xi grid; if a threshold >0 is set, never count ip1=0.
        for(len_t ip1 = iMin; ip1<np1; ip1++){
            for(len_t ip2 = 0; ip2<np2; ip2++){
                const real_t p = mg->GetP(ip1,ip2);
                if(p >= p0){
                    ind = ir*np1*np2+ip2*np1+ip1;
                    this->integrand[ind] = 1;
                }
            }
        }
    }
}

/**
 * Converts the input pThreshold to a threshold in units of m_e*c, depending on which
 * unit was specified in the input (via pThresholdMode).
 */
real_t DensityFromDistributionFunction::GetThreshold(len_t ir, FVM::UnknownQuantityHandler *unknowns){
    switch(pMode){
        case P_THRESHOLD_MODE_MC:
            return pThreshold;
        case P_THRESHOLD_MODE_THERMAL:{
            const real_t Tcold = unknowns->GetUnknownData(OptionConstants::UQTY_T_COLD)[ir];
            return pThreshold * sqrt(2*Tcold/Constants::mc2inEV);
        }
        default:
            throw FVM::FVMException("DensityFromDistributionFunction: Unrecognized p threshold mode.");
            return -1;
    }
}

/**
 * Method that is called whenever the grid is rebuilt. We only
 * need to rebuild this EquationTerm if the total number of grid
 * cells changes.
 */
bool DensityFromDistributionFunction::GridRebuilt() {
    if (this->MomentQuantity::GridRebuilt()) {
        return true;
    } else 
        return false;
}

