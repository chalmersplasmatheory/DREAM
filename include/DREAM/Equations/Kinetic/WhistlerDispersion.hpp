#ifndef _DREAM_EQUATIONS_KINETIC_WHISTLER_DISPERSION_HPP
#define _DREAM_EQUATIONS_KINETIC_WHISTLER_DISPERSION_HPP

#include "DREAM/Constants.hpp"
#include <vector>
#include <string>

namespace DREAM {
    /**
     * Class for solving whistler wave dispersion relations using 
     * PDRF's matrix eigenvalue approach.
     * 
     * This class implements the multi-fluid plasma dispersion relation by constructing
     * a (4S+6)×(4S+6) generalized eigenvalue problem:
     *   M · x = ω · A · x
     * where S is the number of species (typically 2: electrons + ions).
     * 
     * For given k_∥ and k_perp, this directly solves for all wave mode frequencies ω.
     * The approach follows Xie's PDRF code (2014).
     * 
     * Matrix dimensions:
     * - 4 variables per species: (δn, δv_x, δv_y, δv_z)
     * - 6 electromagnetic field variables: (E_x, E_y, E_z, B_x, B_y, B_z)
     * - Total: 4*S + 6 = 14 for S=2 species
     * 
     * Reference:
     * [Xie2014] H. S. Xie, PDRF: A general dispersion relation solver for 
     *            magnetized multi-fluid plasma, CPC 185, 670-675 (2014).
     * Implementation based on: /data/zhzhou/pdrf/pypdrf/pdrf_SI.py
     */
    class WhistlerDispersion {
    private:
        // Plasma parameters
        real_t B0;              // Background magnetic field (Tesla)
        real_t density;         // Electron density (m^-3)
        real_t ion_mass_ratio;  // Ion mass / electron mass (e.g., 1836 for protons)
        real_t Zeff;           // Effective ion charge
        
        // Derived parameters
        real_t omega_ce;       // Electron cyclotron frequency (rad/s)
        real_t omega_ci;       // Ion cyclotron frequency (rad/s)
        real_t omega_pe;       // Electron plasma frequency (rad/s)
        real_t omega_pi;       // Ion plasma frequency (rad/s)
        real_t v_A;            // Alfvén velocity (m/s)
        real_t w_factor;       // Simplified dispersion coefficient: w = ω_ce * c² / ω_pe²
        
        // Method selection
        bool use_simple_dispersion;  // Use simplified whistler dispersion relation
        
        // Physical constants
        static constexpr real_t c = 2.99792458e8;      // Speed of light (m/s)
        static constexpr real_t epsilon0 = 8.854187817e-12;  // Permittivity (F/m)
        static constexpr real_t mu0 = 1.2566370614e-6;     // Permeability (H/m) = 4π×10^-7
        static constexpr real_t e_charge = 1.60217662e-19;  // Elementary charge (C)
        static constexpr real_t m_electron = 9.1094e-31;    // Electron mass (kg)
        
    public:
        /**
         * Constructor
         * 
         * @param B0_val          Background magnetic field (Tesla)
         * @param density_val     Electron density (m^-3)
         * @param ion_mass_ratio_val  Ion to electron mass ratio
         * @param Zeff_val        Effective ion charge
         * @param use_simple      Use simplified whistler dispersion (faster, less accurate)
         */
        WhistlerDispersion(real_t B0_val, real_t density_val, 
                          real_t ion_mass_ratio_val=1836.0, 
                          real_t Zeff_val=1.0,
                          bool use_simple=false);
        
        /**
         * Destructor - cleans up Python resources if this is the last instance
         */
        ~WhistlerDispersion();
        
        /**
         * Calculate wave frequency omega for given (k, theta_k)
         * 
         * Uses PDRF's generalized eigenvalue method to solve the dispersion relation
         * (or simplified formula if use_simple_dispersion is true).
         * For given k and θ, calculates k_x = k*sin(θ), k_z = k*cos(θ),
         * then constructs the (4S+6)×(4S+6) matrix and solves M·x = ω·A·x.
         * Returns the whistler branch solution from all eigenvalues.
         * 
         * @param k       Wavenumber magnitude (m^-1)
         * @param theta   Angle between k and B0 (radians)
         * @return        Wave frequency omega (rad/s), or -1 if no valid solution
         */
        real_t calculateOmega(real_t k, real_t theta) const;
        
        /**
         * Calculate wave frequency using simplified whistler dispersion relation.
         * 
         * For Ω_ci ≪ ω ≪ Ω_ce, the cold plasma whistler dispersion simplifies to:
         *   ω = k|k_| * w,  where w = ω_ce * c² / ω_pe²
         * 
         * This is derived from the cold plasma dielectric tensor in the whistler limit.
         * Reference: calculate_kperp_bon.py (BON code)
         * 
         * @param k       Wavenumber magnitude (m^-1)
         * @param theta   Angle between k and B0 (radians)
         * @return        Wave frequency omega (rad/s)
         */
        real_t calculateOmegaSimple(real_t k, real_t theta) const;
        
        /**
         * Calculate perpendicular wavenumber using simplified whistler dispersion.
         * 
         * Inverse problem: given ω and k_∥, find k_⊥ using:
         *   ω = k|k_∥| * w  =>  k = ω / (|k_∥| * w)
         *   k_⊥ = sqrt(k² - k_∥²)
         * 
         * @param omega   Wave frequency (rad/s)
         * @param k_par   Parallel wavenumber (m^-1)
         * @return        Perpendicular wavenumber k_ (m^-1), or -1 if no real solution
         */
        real_t calculateKperpSimple(real_t omega, real_t k_par) const;
        
        /**
         * Calculate group velocity dω/dk
         * 
         * Uses numerical differentiation for robustness.
         * 
         * @param k       Wavenumber (m^-1)
         * @param theta   Angle (radians)
         * @return        Group velocity (m/s)
         */
        real_t calculateGroupVelocity(real_t k, real_t theta) const;
        
        /**
         * Calculate polarization vectors (Ex, Ey, Ez)
         * 
         * For cold plasma waves, the polarization is determined by the
         * dielectric tensor eigenvectors.
         * 
         * @param omega   Wave frequency (rad/s)
         * @param k       Wavenumber (m^-1)
         * @param theta   Angle (radians)
         * @param Ex_out  Output: Ex component
         * @param Ey_out  Output: Ey component
         * @param Ez_out  Output: Ez component
         */
        void calculatePolarization(real_t omega, real_t k, real_t theta,
                                  real_t &Ex_out, real_t &Ey_out, real_t &Ez_out) const;
        
        /**
         * Calculate the normalization factor for the dispersion relation
         * 
         * This corresponds to the 'den' term in QUADRE, which is related to
         * ∂D/∂ω where D is the dispersion function.
         * 
         * @param omega   Wave frequency (rad/s)
         * @param k       Wavenumber (m^-1)
         * @param theta   Angle (radians)
         * @return        Normalization factor
         */
        real_t calculateDenominator(real_t omega, real_t k, real_t theta) const;
        
        // Accessor methods
        real_t getOmegaCE() const { return omega_ce; }
        real_t getOmegaCI() const { return omega_ci; }
        real_t getOmegaPE() const { return omega_pe; }
        real_t getOmegaPI() const { return omega_pi; }
        real_t getB0() const { return B0; }
        real_t getDensity() const { return density; }
        
        /**
         * Check if frequency is in whistler range
         * 
         * @param omega   Frequency to check (rad/s)
         * @return        True if Ω_ci << ω < Ω_ce
         */
        bool isInWhistlerRange(real_t omega) const;
        
        /**
         * Print dispersion relation information
         */
        void printInfo() const;
        
    private:
        /**
         * Initialize derived parameters from basic inputs
         */
        void initializeParameters();
        
        /**
         * Solve dispersion relation using PDRF's generalized eigenvalue method
         * 
         * Calls Python PDRF solver via subprocess to solve:
         *   M · x = ω · A · x
         * 
         * Returns all valid eigenfrequencies, selects whistler branch.
         * 
         * @param kx      Wave number component perpendicular to B0 (m^-1)
         * @param kz      Wave number component parallel to B0 (m^-1)
         * @return        Vector of valid frequencies (rad/s)
         */
        std::vector<real_t> solvePDRFMethod(real_t kx, real_t kz) const;
        
        /**
         * Solve cold plasma dispersion relation using Stix notation
         * 
         * For cold plasma, the dispersion relation can be written as:
         *   An^4 - Bn^2 + C = 0
         * where n = kc/ω is the refractive index.
         * 
         * Coefficients A, B, C are functions of Stix parameters R, L, P, S, D.
         * 
         * @param k       Wavenumber (m^-1)
         * @param theta   Angle (radians)
         * @return        Refractive index squared n^2
         */
        real_t solveColdPlasmaDispersion(real_t k, real_t theta) const;
        
        /**
         * Calculate Stix parameters (R, L, P, S, D)
         * 
         * @param omega   Wave frequency (rad/s)
         * @param R_out   Right-hand circular polarization
         * @param L_out   Left-hand circular polarization
         * @param P_out   Plasma oscillation term
         * @param S_out   Sum term (R+L)/2
         * @param D_out   Difference term (R-L)/2
         */
        void calculateStixParameters(real_t omega, 
                                    real_t &R_out, real_t &L_out, 
                                    real_t &P_out, real_t &S_out, 
                                    real_t &D_out) const;
    };
}

#endif /* _DREAM_EQUATIONS_KINETIC_WHISTLER_DISPERSION_HPP */
