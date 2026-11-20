#ifndef _DREAM_EQUATIONS_KINETIC_WAVEPARTICLECOUPLING_HPP
#define _DREAM_EQUATIONS_KINETIC_WAVEPARTICLECOUPLING_HPP

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/WhistlerDispersion.hpp"
#include "DREAM/Settings/WaveSpectrum.hpp"
#include <vector>
#include <cmath>

namespace DREAM {
    /**
     * Class for calculating wave-particle coupling strength |Ψ_{n,k}|².
     * 
     * The coupling strength determines the efficiency of energy/momentum
     * transfer between waves and particles at resonance.
     * 
     * Based on quasi-linear theory (from QUADRE implementation):
     *   |Ψ_{n,k}|² = weight / |∂ω/∂k - v_∥|
     *   
     * where weight includes:
     *   - Bessel function terms J_n(k_⊥ρ)
     *   - Wave polarization (Ex, Ey, Ez)
     *   - Normalization by density and plasma parameters
     */
    class WaveParticleCoupling {
    private:
        real_t omega_pe;   // Electron plasma frequency (rad/s)
        real_t omega_ce;   // Electron cyclotron frequency (rad/s)
        real_t omega_pi;   // Ion plasma frequency (rad/s)
        real_t omega_ci;   // Ion cyclotron frequency (rad/s)
        real_t n_e;        // Electron density (m^-3)
        
        // Physical constants
        static constexpr real_t c = 2.99792458e8;      // Speed of light (m/s)
        static constexpr real_t e_charge = 1.60217662e-19;  // Elementary charge (C)
        static constexpr real_t m_electron = 9.1094e-31;    // Electron mass (kg)
        static constexpr real_t m_ion = 3.3436e-27;         // Deuteron mass (kg)
        static constexpr real_t epsilon_0 = 8.854187817e-12; // Vacuum permittivity
        static constexpr real_t delta_omega = 1e-4;         // Small shift for numerical d/dω
        
    public:
        /**
         * Constructor
         * @param omega_pe_val  Electron plasma frequency (rad/s)
         * @param omega_ce_val  Electron cyclotron frequency (rad/s)
         * @param n_e_val       Electron density (m^-3)
         */
        WaveParticleCoupling(real_t omega_pe_val, real_t omega_ce_val, real_t n_e_val);
        
        /**
         * Calculate coupling strength |Ψ_{n,k}|² for given particle state and wave mode.
         * 
         * This implements the formula from QUADRE (inject.py line 689-692):
         *   weight = [n·J_n(k_⊥ρ)/(k_⊥ρ) + E_z·J_n(k_⊥ρ)·ξ/√(1-ξ²) 
         *            - E_y·(J_{n+1}(k_⊥ρ) - J_{n-1}(k_⊥ρ))/2]²
         *   weight = weight / |∂ω/∂k - α|
         *   weight = weight / den
         *   weight = weight · k²·ω_pe² / (n_e·π)
         * 
         * @param p             Particle momentum (normalized to m_e c)
         * @param xi            Pitch-angle cosine (v_∥/v)
         * @param k             Wavenumber (m^-1)
         * @param theta_k       Wave propagation angle (rad)
         * @param n             Harmonic number
         * @param dispersion    Dispersion relation solver
         * @return              Coupling strength |Ψ_{n,k}|² (SI units)
         */
        real_t calculateCouplingStrength(
            real_t p, real_t xi, real_t k, real_t theta_k, int n,
            const WhistlerDispersion &dispersion
        ) const;
        
        /**
         * Calculate bracket term AND dielectric derivative den simultaneously.
         *   bracket = e²/(4π m_e² c²) · |ψ_norm|² / |v_g - v∥ cos θ_k|
         *   den = ∂D/∂ω · ω  (from cold plasma dielectric tensor)
         *
         * @param bracket   [out] The normalized bracket term
         * @param den       [out] Dielectric derivative ∂D/∂ω · ω
         */
        void calculateBracketAndDen(
            real_t p, real_t xi, real_t k, real_t theta_k, int n,
            const WhistlerDispersion &dispersion,
            real_t &bracket, real_t &den
        ) const;
        
        /**
         * Calculate wave polarization (Ex, Ey, Ez) and dielectric derivative (den).
         * Based on cold plasma dielectric tensor for whistler wave.
         *
         * @param omega     Wave angular frequency (rad/s)
         * @param k         Wavenumber (m^-1)
         * @param theta_k   Wave propagation angle (rad)
         * @param Ex        [out] x-component of wave electric field (normalized)
         * @param Ey        [out] y-component of wave electric field
         * @param Ez        [out] z-component of wave electric field
         * @param den       [out] ∂D/∂ω · ω, where D is the dielectric dispersion function
         */
        void calculate_den(
            real_t omega, real_t k, real_t theta_k,
            real_t &Ex, real_t &Ey, real_t &Ez, real_t &den
        ) const;
        
    private:
        /**
         * Calculate Bessel function J_n(x).
         * Uses GSL or standard library implementation.
         */
        real_t besselJ(int n, real_t x) const;
        
        /**
         * Calculate group velocity ∂ω/∂k.
         * Uses numerical differentiation of the dispersion relation.
         */
        real_t calculateGroupVelocity(
            real_t k, real_t theta_k,
            const WhistlerDispersion &dispersion,
            real_t dk = 0.01
        ) const;
        
        /**
         * Calculate gyroradius ρ = p_⊥/(m_e Ω_ce).
         */
        real_t calculateGyroradius(real_t p, real_t xi) const;
    };
}

#endif /* _DREAM_EQUATIONS_KINETIC_WAVEPARTICLECOUPLING_HPP */
