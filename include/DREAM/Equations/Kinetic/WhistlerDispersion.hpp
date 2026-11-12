#ifndef _DREAM_EQUATIONS_KINETIC_WHISTLER_DISPERSION_HPP
#define _DREAM_EQUATIONS_KINETIC_WHISTLER_DISPERSION_HPP

#include "DREAM/Constants.hpp"
#include <vector>
#include <string>

namespace DREAM {
    /**
     * Class for solving whistler wave dispersion relations using 
     * the full cold plasma Stix dielectric tensor.
     * 
     * Solves the Appleton-Hartree dispersion relation:
     *   A·n⁴ - B·n² + C = 0
     * where n = kc/ω is the refractive index and A,B,C depend on
     * Stix parameters S, D, P.
     * 
     * For given (k, θ), uses numerical root-finding to invert the
     * dispersion relation and find ω.
     */
    class WhistlerDispersion {
    private:
        // Plasma parameters
        real_t B0;              // Background magnetic field (Tesla)
        real_t density;         // Electron density (m^-3)
        real_t ion_mass_ratio;  // Ion mass / electron mass (e.g., 1836 for protons)
        real_t Zeff;            // Effective ion charge
        
        // Derived parameters
        real_t omega_ce;       // Electron cyclotron frequency (rad/s)
        real_t omega_ci;       // Ion cyclotron frequency (rad/s)
        real_t omega_pe;       // Electron plasma frequency (rad/s)
        real_t omega_pi;       // Ion plasma frequency (rad/s)
        real_t v_A;            // Alfvén velocity (m/s)
        real_t w_factor;       // Simplified dispersion coefficient (legacy)
        
        // Method selection
        bool use_simple_dispersion;
        
        // Physical constants
        static constexpr real_t c = 2.99792458e8;
        static constexpr real_t epsilon0 = 8.854187817e-12;
        static constexpr real_t mu0 = 1.2566370614e-6;
        static constexpr real_t e_charge = 1.60217662e-19;
        static constexpr real_t m_electron = 9.1094e-31;
        
    public:
        WhistlerDispersion(real_t B0_val, real_t density_val, 
                          real_t ion_mass_ratio_val=1836.0, 
                          real_t Zeff_val=1.0,
                          bool use_simple=false);
        
        ~WhistlerDispersion();
        
        /**
         * Calculate wave frequency omega for given (k, theta_k)
         * Uses full Stix tensor with numerical inversion, or simplified
         * formula if use_simple_dispersion=true.
         */
        real_t calculateOmega(real_t k, real_t theta) const;
        
        /**
         * Simplified whistler dispersion: ω = k|k_∥| * w
         */
        real_t calculateOmegaSimple(real_t k, real_t theta) const;
        
        /**
         * Inverse simplified dispersion: given ω and k_∥, find k_⊥
         */
        real_t calculateKperpSimple(real_t omega, real_t k_par) const;
        
        real_t calculateGroupVelocity(real_t k, real_t theta) const;
        
        void calculatePolarization(real_t omega, real_t k, real_t theta,
                                  real_t &Ex_out, real_t &Ey_out, real_t &Ez_out) const;
        
        real_t calculateDenominator(real_t omega, real_t k, real_t theta) const;
        
        // Accessors
        real_t getOmegaCE() const { return omega_ce; }
        real_t getOmegaCI() const { return omega_ci; }
        real_t getOmegaPE() const { return omega_pe; }
        real_t getOmegaPI() const { return omega_pi; }
        real_t getB0() const { return B0; }
        real_t getDensity() const { return density; }
        
        bool isInWhistlerRange(real_t omega) const;
        void printInfo() const;
        
    private:
        void initializeParameters();
        
        /**
         * Calculate Stix parameters S, D, P for cold plasma
         */
        void calculateStixParameters(real_t omega, 
                                    real_t &S, real_t &D, real_t &P) const;
        
        /**
         * Solve Appleton-Hartree quadratic A·n⁴ - B·n² + C = 0
         * Returns all valid positive real roots for n²
         */
        std::vector<real_t> solveRefractiveIndex(real_t omega, real_t theta) const;
    };
}

#endif /* _DREAM_EQUATIONS_KINETIC_WHISTLER_DISPERSION_HPP */