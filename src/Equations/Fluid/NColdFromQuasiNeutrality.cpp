/**
 * Calculate the cold electron density from quasi-neutrality, i.e.
 * by summing over all ion densities, multiplying with their charges.
 */

#include <iostream>
#include "DREAM/Equations/Fluid/NColdFromQuasiNeutrality.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
NColdFromQuasiNeutrality::NColdFromQuasiNeutrality(
    FVM::Grid *grid, IonHandler *ihdl,
    const len_t idNHot, const len_t idNRE
) : PredeterminedParameter(grid), ions(ihdl), id_nhot(idNHot), id_nre(idNRE) { }

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
    const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    ions->evaluateFreeElectronDensityFromQuasiNeutrality(this->currentData);

    // Subtract n_hot and n_re
    const real_t *nhot = unknowns->GetUnknownData(this->id_nhot);
    const real_t *nre  = unknowns->GetUnknownData(this->id_nre);
    const len_t nr     = this->grid->GetNr();

    for (len_t i = 0; i < nr; i++)
        this->currentData[i] -= nhot[i];
    for (len_t i = 0; i < nr; i++)
        this->currentData[i] -= nre[i];
}

