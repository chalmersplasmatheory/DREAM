#ifndef _DREAM_EQUATIONS_MOLLER_PRODUCTION_RATE_TERM_HPP
#define _DREAM_EQUATIONS_MOLLER_PRODUCTION_RATE_TERM_HPP

#include <cmath>
#include <limits>
#include <vector>

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/config.h"

namespace DREAM {

/**
 * Fluid equation term for the total Møller knock-on production rate.
 *
 * Represents a scalar source proportional to n_tot(r):
 *   dn_re/dt  ~  n_tot(r) * ∫ dp1 dxi1 Vp1 f_primary(p1,xi1) * sigmaTot(p1)
 *
 * This term is explicit in time in f_primary: it reads f_primary at the previous
 * time level in Rebuild() and caches a per-radius rate[] which is used by
 * SetVectorElements() and SetJacobianBlock().
 *
 * The supplied MollerEnergyKernel is used as the single source of truth for
 * sigmaTot(p1) (TotalCS(k)).
 */
class MollerProductionRateTerm : public FVM::EquationTerm {
   private:
    const FVM::Grid *gridTarget = nullptr;   // where secondaries land (for sigmaTot definition)
    const FVM::Grid *gridPrimary = nullptr;  // where primaries live
    FVM::UnknownQuantityHandler *unknowns = nullptr;

    len_t id_f_primary;
    len_t id_ntot;
    const MollerEnergyKernel *energyKernel = nullptr;

    real_t scaleFactor = 1.0;

    // Cached per-radius production rate (without n_tot factor): size Nr(fluid grid)
    std::vector<real_t> rate;
    real_t t_rebuilt = -std::numeric_limits<real_t>::infinity();

    void Validate() const;

   public:
    MollerProductionRateTerm(
        FVM::Grid *fluidGrid, const FVM::Grid *grid_target, const FVM::Grid *grid_primary,
        FVM::UnknownQuantityHandler *u, len_t id_f_primary, const MollerEnergyKernel *energy_kernel,
        real_t scaleFactor = 1.0
    );

    virtual ~MollerProductionRateTerm() = default;

    virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *uqh)
        override;
    virtual void SetVectorElements(real_t *vec, const real_t *n_tot) override;

    virtual bool SetJacobianBlock(
        const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x
    ) override;

    virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
};

}  // namespace DREAM

#endif
