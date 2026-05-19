#ifndef _DREAM_EQUATIONS_KINETIC_QUASILINEARDIFFUSIONTERM_HPP
#define _DREAM_EQUATIONS_KINETIC_QUASILINEARDIFFUSIONTERM_HPP

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/ResonanceSolver.hpp"
#include "DREAM/Equations/Kinetic/WaveParticleCoupling.hpp"
#include "DREAM/Equations/Kinetic/WhistlerDispersion.hpp"
#include "DREAM/Settings/WaveSpectrum.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include <vector>

// Forward declaration for pre-computed matrix loader
namespace DREAM {
namespace Equations {
namespace Kinetic {
class QLMatrixLoader;
}
}
}

namespace DREAM {
    /**
     * Quasi-linear diffusion term from wave-particle interactions.
     * 
     * Implements the diffusion operator:
     *   C[f] = ∂/∂p (D_pp ∂f/∂p + D_pξ ∂f/∂ξ)
     *        + 1/p² ∂/∂ξ (D_pξ p² ∂f/∂p + D_ξξ p² ∂f/∂ξ)
     * 
     * where diffusion coefficients are:
     *   D_pp  = Σ_{n,k} |Ψ_{n,k}|²
     *   D_pξ  = Σ_{n,k} |Ψ_{n,k}|² · (k_∥ v_⊥ / ω)
     *   D_ξξ  = Σ_{n,k} |Ψ_{n,k}|² · (k_∥ v_⊥ / ω)²
     * 
     * The sum is over all wave modes and harmonic numbers at resonance.
     */
    class QuasilinearDiffusionTerm : public FVM::DiffusionTerm {
    private:
        // Traditional components (for on-the-fly computation)
        WaveSpectrum *spectrum;
        WhistlerDispersion *dispersion;
        ResonanceSolver *resonanceSolver;
        WaveParticleCoupling *coupling;
        
        std::vector<int> harmonicModes;  // List of n values
        
        // Pre-computed matrix loader (alternative to on-the-fly computation)
        Equations::Kinetic::QLMatrixLoader *matrix_loader;
        bool use_precomputed_matrix;
        
        // Cached dispersion relation results (omega for each mode)
        real_t *omega_cache;  // [numModes]
        bool dispersion_cached;
        
        // Cached diffusion coefficients (for performance)
        real_t *D_pp_cache;   // [nr][np1+1][np2]
        real_t *D_pxi_cache;  // [nr][np1][np2+1]
        real_t *D_xixi_cache; // [nr][np1][np2+1]
        
        len_t nr_cached, np1_cached, np2_cached;
        bool first_rebuild_done;  // Flag: true after first Rebuild() call
        bool operator_cached;     // Flag: true if diffusion operator matrices are cached
        
        // Track wave amplitudes to detect changes
        real_t *last_amplitudes;  // [numModes] - stores amplitudes from last calculation
        len_t num_modes_tracked;
        
        // Time-dependent amplitude (for pre-computed matrix scaling)
        real_t current_amplitude;
        
        /**
         * Quick resonance check: determine if a grid point can possibly resonate
         * with any wave mode. This is a fast pre-filter before calling the full
         * resonance solver.
         * 
         * Based on QUADRE's approach:
         *   omega_res = k*v_parallel*cos(theta_k) + n*omega_ce/gamma
         *   
         * For a given (p, xi), calculate the range of possible resonant frequencies
         * as k varies over [k_min, k_max]. If this range doesn't overlap with
         * the wave spectrum frequency range, skip this point.
         * 
         * @param p          Momentum
         * @param xi         Pitch-angle cosine
         * @param theta_k    Wave propagation angle
         * @param n          Harmonic number
         * @return           true if resonance is possible, false otherwise
         */
        bool CanResonate(
            real_t p,
            real_t xi,
            real_t theta_k,
            int n
        ) const;
        
    public:
        /**
         * Constructor (traditional on-the-fly computation)
         */
        QuasilinearDiffusionTerm(
            FVM::Grid *grid,
            WaveSpectrum *spectrum,
            WhistlerDispersion *dispersion,
            ResonanceSolver *resonanceSolver,
            WaveParticleCoupling *coupling,
            const std::vector<int> &harmonicModes
        );
        
        /**
         * Constructor (using pre-computed matrix from HDF5)
         * @param grid              Momentum grid
         * @param hdf5_file         Path to pre-computed HDF5 file
         * @param initial_amplitude Initial wave amplitude (can be updated at runtime)
         */
        QuasilinearDiffusionTerm(
            FVM::Grid *grid,
            const std::string &hdf5_file,
            real_t initial_amplitude = 1.0
        );
        
        ~QuasilinearDiffusionTerm();
        
        /**
         * Rebuild diffusion coefficients
         */
        void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *uqh) override;
        
        /**
         * Set current wave amplitude (for pre-computed matrix mode)
         * @param A_t  Current amplitude at time t
         */
        void setCurrentAmplitude(real_t A_t);
        
        /**
         * Get current wave amplitude
         * @return Current amplitude A(t)
         */
        real_t getCurrentAmplitude() const { return current_amplitude; }
        
    private:
        /**
         * Pre-calculate and cache dispersion relation for all modes
         */
        void precalculateDispersion();
        
        /**
         * Calculate diffusion coefficients for given radial grid point
         */
        void calculateDiffusionCoefficients(len_t ir);
        
        // New method using reverse resonance solving
        void calculateDiffusionCoefficientsReverse(len_t ir);
        
        /**
         * Get cached diffusion coefficient D_pp at grid point
         */
        real_t getD_pp(len_t ir, len_t i, len_t j) const;
        
        /**
         * Get cached diffusion coefficient D_pξ at grid point
         */
        real_t getD_pxi(len_t ir, len_t i, len_t j) const;
        
        /**
         * Get cached diffusion coefficient D_ξξ at grid point
         */
        real_t getD_xixi(len_t ir, len_t i, len_t j) const;
    };
}

#endif /* _DREAM_EQUATIONS_KINETIC_QUASILINEARDIFFUSIONTERM_HPP */
