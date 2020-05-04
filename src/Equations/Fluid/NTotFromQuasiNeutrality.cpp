/**
 * Calculate the total electron density from quasi-neutrality, i.e.
 * by summing over all ion densities, multiplying with their charges.
 */

#include "DREAM/Equations/Fluid/NTotFromQuasiNeutrality.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
NTotFromQuasiNeutrality::NTotFromQuasiNeutrality(
    FVM::Grid *grid, IonHandler *ihdl
) : PredeterminedParameter(grid), ions(ihdl) { }

/**
 * Destructor.
 */
NTotFromQuasiNeutrality::~NTotFromQuasiNeutrality() { }


/**
 * Recalculate the total electron density.
 *
 * (All input parameters are unused)
 * t:        Current (simulation) time.
 * dt:       Length of time step to take.
 * unknowns: List of unknown quantities being solved for.
 */
void NTotFromQuasiNeutrality::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler*
) {
    ions->evaluateFreePlusBoundElectronDensityFromQuasiNeutrality(this->currentData);
}

