#ifndef _DREAM_EQUATIONS_KINETIC_RESONANCE_SOLVER_HPP
#define _DREAM_EQUATIONS_KINETIC_RESONANCE_SOLVER_HPP

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/WhistlerDispersion.hpp"
#include <vector>
#include <cmath>

namespace DREAM {
    /**
     * Class for solving the resonance condition between waves and particles.
     * 
     * The resonance condition is:
     *   ω - k_∥ v_∥ + n Ω_ce / γ = 0
     * 
     * Combined with the dispersion relation ω = ω(k, θ_k), this becomes
     * a nonlinear equation for k that must be solved numerically.
     * 
     * This class uses root-finding methods (brentq) to find resonant
     * wavenumbers for given particle states and wave angles.
     */
    class ResonanceSolver {
    private:
        real_t omega_ce;   // Electron cyclotron frequency (rad/s)
        
        // Physical constants
        static constexpr real_t c = 2.99792458e8;      // Speed of light (m/s)
        static constexpr real_t e_charge = 1.60217662e-19;  // Elementary charge (C)
        static constexpr real_t m_electron = 9.1094e-31;    // Electron mass (kg)
        
    public:
        /**
         * Constructor
         * @param omega_ce_val  Electron cyclotron frequency (rad/s)
         */
        ResonanceSolver(real_t omega_ce_val);
        
        /**
         * Find resonant wavenumber(s) for given particle state and wave angle.
         * 
         * Solves: ω(k, θ_k) = k cos(θ_k) v_∥ - n Ω_ce / γ
         * 
         * Uses numerical root-finding to locate k values where the resonance
         * condition is satisfied. May return multiple solutions if they exist.
         * 
         * @param p             Particle momentum (normalized to m_e c)
         * @param xi            Pitch-angle cosine (v_∥ / v)
         * @param theta_k       Wave propagation angle (rad)
         * @param n             Harmonic number (typically -1, 0, or +1)
         * @param dispersion    Dispersion relation solver
         * @param k_min         Minimum search wavenumber (m^-1)
         * @param k_max         Maximum search wavenumber (m^-1)
         * @return              Vector of resonant wavenumbers (may be empty)
         */
        std::vector<real_t> findResonantK(
            real_t p, real_t xi, real_t theta_k, int n,
            const WhistlerDispersion &dispersion,
            real_t k_min = 1.0, real_t k_max = 200.0
        ) const;
        
        /**
         * Find resonant wavenumber(s) using pre-calculated omega values.
         * This version avoids repeated PDRF calls by using cached dispersion data.
         * 
         * @param p             Particle momentum (normalized to m_e c)
         * @param xi            Pitch-angle cosine (v_∥ / v)
         * @param theta_k       Wave propagation angle (rad)
         * @param n             Harmonic number
         * @param omega_cache   Pre-calculated omega values for each mode
         * @param k_values      Corresponding k values for each mode
         * @param num_modes     Number of modes in the cache
         * @return              Vector of resonant wavenumbers (may be empty)
         */
        std::vector<real_t> findResonantKWithCachedOmega(
            real_t p, real_t xi, real_t theta_k, int n,
            const real_t* omega_cache,
            const real_t* k_values,
            len_t num_modes
        ) const;
        
        /**
         * Check if a particle-wave pair is in resonance.
         * 
         * @param omega         Wave frequency (rad/s)
         * @param kpar          Parallel wavenumber k_∥ (m^-1)
         * @param vpar          Parallel velocity (m/s)
         * @param gamma         Lorentz factor
         * @param n             Harmonic number
         * @param tolerance     Relative tolerance for resonance check
         * @return              True if |ω - k_∥ v_∥ + n Ω_ce/γ| / ω < tolerance
         */
        bool isInResonance(
            real_t omega, real_t kpar, real_t vpar,
            real_t gamma, int n, real_t tolerance = 1e-3
        ) const;
        
        /**
         * Find resonant momentum p for given wave parameters and pitch angle.
         * This is the REVERSE approach: instead of checking if a grid point is resonant,
         * we solve for the exact resonant p given (k, theta_k, n, xi).
         * 
         * Solves: ω(k, θ_k) - k cos(θ_k) v_∥(p,xi) + n Ω_ce / γ(p) = 0
         * 
         * For whistler waves with simplified dispersion ω = k|k_∥| * w:
         *   k|k_∥| * w - k cos(θ_k) * (p/γ) * xi * c + n Ω_ce / γ = 0
         * 
         * This can be rearranged to a quadratic equation in p.
         * 
         * @param k             Wavenumber (m^-1)
         * @param theta_k       Wave propagation angle (rad)
         * @param n             Harmonic number
         * @param xi            Pitch-angle cosine
         * @param dispersion    Dispersion relation solver
         * @param p_min         Minimum search momentum
         * @param p_max         Maximum search momentum
         * @return              Vector of resonant momenta (may be empty)
         */
        std::vector<real_t> findResonantP(
            real_t k, real_t theta_k, int n, real_t xi,
            const WhistlerDispersion &dispersion,
            real_t p_min = 0.0, real_t p_max = 100.0
        ) const;
        
    private:
        /**
         * Calculate parallel velocity from momentum and pitch angle.
         * @param p     Momentum (normalized to m_e c)
         * @param xi    Pitch-angle cosine
         * @return      v_∥ / c
         */
        real_t calculateVparallel(real_t p, real_t xi) const;
        
        /**
         * Calculate Lorentz factor from momentum.
         * @param p     Momentum (normalized to m_e c)
         * @return      γ = sqrt(1 + p^2)
         */
        real_t calculateGamma(real_t p) const;
        
        /**
         * Resonance residual function.
         * Returns: ω_disp(k, θ_k) - [k cos(θ_k) v_∥ - n Ω_ce / γ]
         * 
         * Root of this function corresponds to resonance condition.
         */
        real_t resonanceResidual(
            real_t k, real_t vpar_over_c, real_t gamma,
            real_t cos_theta, int n,
            const WhistlerDispersion &dispersion, real_t theta_k
        ) const;
    };
}

#endif /* _DREAM_EQUATIONS_KINETIC_RESONANCE_SOLVER_HPP */
