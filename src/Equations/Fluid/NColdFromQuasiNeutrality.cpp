/**
 * Calculate the cold electron density from quasi-neutrality, i.e.
 * by summing over all ion densities, multiplying with their charges.
 */

#include "DREAM/Equations/Fluid/NColdFromQuasiNeutrality.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
NColdFromQuasiNeutrality::NColdFromQuasiNeutrality(
    FVM::Grid *grid, IonHandler *ihdl
) : PredeterminedParameter(grid), ions(ihdl) { }

/**
 * Destructor.
 */
NColdFromQuasiNeutrality::~NColdFromQuasiNeutrality() { }


/**
 * Recalculate the cold electron density.
 *
 * (All input parameters are unused)
 * t:        Current (simulation) time.
 * dt:       Length of time step to take.
 * unknowns: List of unknown quantities being solved for.
 */
void NColdFromQuasiNeutrality::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler*
) {
    ions->evaluateFreeElectronDensityFromQuasiNeutrality(this->currentData);
}

