#ifndef _DREAM_EQUATIONS_MOLLER_ENERGY_KERNEL_HPP
#define _DREAM_EQUATIONS_MOLLER_ENERGY_KERNEL_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/config.h"

namespace DREAM {

/**
 * Helper which owns the (i,k) Møller momentum-space kernel integrated
 * over outgoing p-cells:
 *   S_{ik} = ∫_{cell i} dσ/dp ...   (as provided by KnockOnUtilities)
 */
class MollerEnergyKernel {
   private:
    const FVM::Grid *gridK = nullptr;
    const FVM::Grid *gridP = nullptr;
    real_t pCutoff = 0;

    // table of cell-averaged differential cross sections
    real_t *Sik = nullptr;  // size: NpK0 * NpP0
    len_t NpK0 = 0;
    len_t NpP0 = 0;

    void Deallocate();
    void ValidateInputParameters() const;
    void ValidateGridAssumptions() const;

   public:
    MollerEnergyKernel(
        const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t p_cutoff
    );
    ~MollerEnergyKernel();

    void GridRebuilt();

    // Cell-integrated (over knock-on momenta) Møller differential cross section weighted by v1.
    // Momentum cell i (clamped by pCutoff and half-max energy transfer) on the
    // knock-on grid, and p1_k on primary grid.
    real_t DifferentialCS(len_t i, len_t k) const {
        // Same across radii under validated assumptions.
        return Sik[i * NpP0 + k];
    }
    real_t *DifferentialCS() const {
        // Same across radii under validated assumptions.
        return Sik;
    }

    // The total cross section for collisions between the knock-on grid above pCutoff,
    // evaluated at momentum p1_k on the primary grid.
    real_t TotalCS(len_t k) const;
};

}  // namespace DREAM

#endif
