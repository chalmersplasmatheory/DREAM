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
     *   ω - k_∥ v_∥ - n Ω_ce / γ = 0
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
         * Find resonant momentum p for given wave parameters and pitch angle.
         * This is the REVERSE approach: instead of checking if a grid point is resonant,
         * we solve for the exact resonant p given (k, theta_k, n, xi).
         * 
         * Solves: ω(k, θ_k) - k cos(θ_k) v_∥(p,xi) - n Ω_ce / γ(p) = 0
         * 
         * For whistler waves with simplified dispersion ω = k|k_∥| * w:
         *   k|k_∥| * w - k cos(θ_k) * (p/γ) * xi * c - n Ω_ce / γ = 0
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
        // Note: Helper functions calculateVparallel, calculateGamma, and resonanceResidual
        // were removed as they were only used by deprecated methods (findResonantK, etc.)
    };
}

#endif /* _DREAM_EQUATIONS_KINETIC_RESONANCE_SOLVER_HPP */
