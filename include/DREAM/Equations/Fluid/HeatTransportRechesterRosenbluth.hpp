#ifndef _DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RECHESTER_ROSENBLUTH_HPP
#define _DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RECHESTER_ROSENBLUTH_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

// Collision-related quantities are accessed via RunawayFluid
#include "DREAM/Equations/RunawayFluid.hpp"

namespace DREAM {

    /**
     * Rechester–Rosenbluth heat transport operator with
     * collision-modified effective diffusivity.
     *
     * This class implements a generalized RR-type heat diffusion model,
     * including:
     *  - stochastic magnetic field transport (via deltaB/B),
     *  - collision factor eta,
     *  - multiple competing transport channels combined via a soft-min.
     */
    class HeatTransportRechesterRosenbluth : public FVM::DiffusionTerm {
    protected:
        // Interpolator providing δB/B on the radial (flux) grid
        FVM::Interpolator1D *deltaBOverB = nullptr;

        // Handler for accessing unknown quantities (n, T, j, ...)
        FVM::UnknownQuantityHandler *unknowns = nullptr;

        // Base RR diffusion prefactor D (defined on flux grid, size nr+1)
        real_t *dD = nullptr;

        // Cache of δB/B values on flux grid, stored for Jacobian evaluation
        real_t *dBB_cache = nullptr;

        /**
         * Soft-min hardness parameter.
         * Larger pSoft → sharper transition (closer to hard minimum)
         * Smaller pSoft → smoother interpolation between regimes
         */
        real_t pSoft = 4.0;

        // ------------------------------------------------------------------
        // Caches for competing transport coefficients (defined on flux grid)
        // ------------------------------------------------------------------

        // chiQL: collision-modified quasi-linear RR diffusivity
        real_t *chiQL_cache  = nullptr;

        // chi2: strong-turbulence / short-connection-length asymptote
        real_t *chi2_cache   = nullptr;

        // chi3: collision-limited (ballistic / finite-mean-free-path) channel
        real_t *chi3_cache   = nullptr;

        // chiEff: effective diffusivity obtained from soft-min combination
        real_t *chiEff_cache = nullptr;

        // ------------------------------------------------------------------
        // IDs of unknown quantities used by this operator
        // ------------------------------------------------------------------

        len_t id_n_cold;   // cold electron density
        len_t id_T_cold;   // cold electron temperature
        len_t id_j_tot;    // total plasma current (used for q-profile)

        // ------------------------------------------------------------------
        // Additional members for collision-modified eta model
        // ------------------------------------------------------------------

        // Access to thermal collision time τ_th via RunawayFluid
        RunawayFluid *REFluid = nullptr;

        // Parallel correlation length scale (typically ~ qR0)
        real_t L0  = 0.0;

        // Critical perpendicular correlation length (used for L_perp)
        real_t Lcr = 0.0;

        // --- Caches for eta model (defined on flux grid) ---

        // eta: multiplicative suppression factor for RR transport
        real_t *eta_cache     = nullptr;

        // LK: Kolmogorov (or Kubo-like) connection length
        real_t *LK_cache      = nullptr;

        // LKdelta: logarithmically corrected LK entering eta
        real_t *LKdelta_cache = nullptr;

        // ------------------------------------------------------------------

        /**
         * Allocate and initialize diffusion coefficient arrays and caches.
         * Called at construction and whenever the grid is rebuilt.
         */
        void AllocateDiffCoeff();

        /**
         * Set partial derivatives of the diffusion coefficient
         * with respect to unknown quantities (Jacobian contribution).
         */
        virtual void SetPartialDiffusionTerm(len_t, len_t) override;

        /**
         * Evaluate δB/B at time t using the provided interpolator.
         */
        virtual const real_t *EvaluateDeltaBOverB(const real_t t) {
            return this->deltaBOverB->Eval(t);
        }

    public:
        /**
         * Constructor.
         *
         * @param grid        Computational grid
         * @param dB_B        Interpolator for δB/B
         * @param unknowns   Unknown quantity handler
         * @param REFluid    Access to collision times via runaway-fluid module
         * @param L0         Parallel correlation length scale
         * @param Lcr        Critical perpendicular correlation length
         */
        HeatTransportRechesterRosenbluth(
            FVM::Grid*,
            FVM::Interpolator1D*,
            FVM::UnknownQuantityHandler*,
            RunawayFluid* REFluid,
            real_t L0,
            real_t Lcr
        );

        /**
         * Destructor.
         */
        virtual ~HeatTransportRechesterRosenbluth();

        /**
         * Called when the grid is rebuilt.
         * Reallocates all coefficient arrays on the new grid.
         */
        virtual bool GridRebuilt() override;

        /**
         * Recompute diffusion coefficients at the current time.
         */
        virtual void Rebuild(
            const real_t,
            const real_t,
            FVM::UnknownQuantityHandler*
        ) override;
    };

}

#endif /* _DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RECHESTER_ROSENBLUTH_HPP */
