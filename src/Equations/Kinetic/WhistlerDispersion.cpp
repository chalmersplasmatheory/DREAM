/**
 * WhistlerDispersion class - Full Stix Dielectric Tensor Implementation
 * 
 * Solves the exact cold plasma dispersion relation using Stix parameters (S, D, P)
 * and the Appleton-Hartree quadratic equation for the refractive index.
 * 
 */

#include "DREAM/Equations/Kinetic/WhistlerDispersion.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

using namespace DREAM;

/**
 * Constructor
 */
WhistlerDispersion::WhistlerDispersion(real_t B0_val, real_t density_val, 
                                      real_t ion_mass_ratio_val, 
                                      real_t Zeff_val,
                                      bool use_simple)
    : B0(B0_val), density(density_val), 
      ion_mass_ratio(ion_mass_ratio_val), Zeff(Zeff_val),
      use_simple_dispersion(use_simple) {
    
    initializeParameters();
}

/**
 * Destructor
 */
WhistlerDispersion::~WhistlerDispersion() {
}

/**
 * Initialize derived physical parameters
 */
void WhistlerDispersion::initializeParameters() {
    // Cyclotron frequencies (positive magnitudes)
    omega_ce = e_charge * B0 / m_electron;
    omega_ci = omega_ce / ion_mass_ratio;
    
    // Plasma frequencies
    omega_pe = std::sqrt(density * e_charge * e_charge / (epsilon0 * m_electron));
    
    // Ion density from quasi-neutrality: n_e = Zeff * n_i
    real_t n_ion = density / Zeff;
    real_t m_ion = ion_mass_ratio * m_electron;
    // Assuming singly charged ions for plasma frequency (Z=1 for density calc, Zeff for collisions)
    // If ions have charge Z_i, omega_pi^2 = n_i * (Z_i * e)^2 / (epsilon0 * m_i)
    // For simplicity and standard Stix, we assume Z_i = 1 for the ion species contributing to dispersion
    omega_pi = std::sqrt(n_ion * e_charge * e_charge / (epsilon0 * m_ion));
    
    // Alfvén velocity
    v_A = B0 / std::sqrt(mu0 * n_ion * m_ion);
    
    // Legacy w_factor (kept for backward compatibility if use_simple_dispersion=true)
    w_factor = omega_ce * c * c / (omega_pe * omega_pe);
    
    if (use_simple_dispersion) {
        std::cout << "WhistlerDispersion initialized (SIMPLIFIED dispersion relation):" << std::endl;
        std::cout << "  Dispersion: ω = k|k_∥| * w" << std::endl;
    } else {
        std::cout << "WhistlerDispersion initialized (Full Stix Dielectric Tensor):" << std::endl;
        std::cout << "  B0 = " << B0 << " T" << std::endl;
        std::cout << "  n_e = " << density << " m^-3" << std::endl;
        std::cout << "  ω_ce = " << omega_ce/(2*M_PI)/1e9 << " GHz" << std::endl;
        std::cout << "  ω_pe = " << omega_pe/(2*M_PI)/1e9 << " GHz" << std::endl;
        std::cout << "  ω_ci = " << omega_ci/(2*M_PI)/1e6 << " MHz" << std::endl;
    }
}

/**
 * Calculate Stix parameters S, D, P for a given frequency omega.
 * Cold plasma approximation.
 */
void WhistlerDispersion::calculateStixParameters(real_t omega, real_t &S, real_t &D, real_t &P) const {
    real_t omega2 = omega * omega;
    
    // Electron contributions
    real_t denom_e = omega2 - omega_ce * omega_ce;
    // Avoid exact resonance singularity (though whistler is below omega_ce)
    if (std::abs(denom_e) < 1e-10 * omega2) denom_e = 1e-10 * omega2; 
    
    // Ion contributions
    real_t denom_i = omega2 - omega_ci * omega_ci;
    if (std::abs(denom_i) < 1e-10 * omega2) denom_i = 1e-10 * omega2;

    // S = 1 - sum( omega_ps^2 / (omega^2 - Omega_s^2) )
    S = 1.0 - (omega_pe * omega_pe) / denom_e - (omega_pi * omega_pi) / denom_i;
    
    // D = sum( (Omega_s / omega) * omega_ps^2 / (omega^2 - Omega_s^2) )
    // Note: Stix convention with positive cyclotron frequencies requires careful sign handling.
    // For electrons (charge -e), the contribution to D is negative.
    // For ions (charge +e), the contribution to D is positive.
    D = - (omega_ce / omega) * (omega_pe * omega_pe) / denom_e 
        + (omega_ci / omega) * (omega_pi * omega_pi) / denom_i;
        
    // P = 1 - sum( omega_ps^2 / omega^2 )
    P = 1.0 - (omega_pe * omega_pe) / omega2 - (omega_pi * omega_pi) / omega2;
}

/**
 * Calculate the refractive index squared (n^2) for a given omega and propagation angle theta.
 * Solves the Appleton-Hartree quadratic equation: A * n^4 - B * n^2 + C = 0
 * 
 * Returns the two roots. For Whistler waves, we typically want the fast/whistler branch.
 */
std::vector<real_t> WhistlerDispersion::solveRefractiveIndex(real_t omega, real_t theta) const {
    std::vector<real_t> n2_roots;
    
    real_t S, D, P;
    calculateStixParameters(omega, S, D, P);
    
    real_t sin_theta = std::sin(theta);
    real_t cos_theta = std::cos(theta);
    real_t sin2 = sin_theta * sin_theta;
    real_t cos2 = cos_theta * cos_theta;
    
    // Coefficients of the quadratic equation A*n^4 - B*n^2 + C = 0
    real_t A = S * sin2 + P * cos2;
    real_t B = (S * S - D * D) * sin2 + P * S * (1.0 + cos2);
    real_t C = P * (S * S - D * D);
    
    // Handle degenerate cases
    if (std::abs(A) < 1e-12) {
        if (std::abs(B) > 1e-12) {
            n2_roots.push_back(C / B);
        }
        return n2_roots;
    }
    
    real_t discriminant = B * B - 4.0 * A * C;
    
    if (discriminant < 0) {
        // No real roots (evanescent wave)
        return n2_roots;
    }
    
    real_t sqrt_disc = std::sqrt(discriminant);
    
    // Two roots for n^2
    real_t n2_1 = (B + sqrt_disc) / (2.0 * A);
    real_t n2_2 = (B - sqrt_disc) / (2.0 * A);
    
    if (n2_1 > 0) n2_roots.push_back(n2_1);
    if (n2_2 > 0 && std::abs(n2_2 - n2_1) > 1e-10) n2_roots.push_back(n2_2);
    
    return n2_roots;
}

/**
 * Calculate wave frequency omega for given (k, theta) using the full Stix tensor.
 * 
 * This is an inverse problem: given k and theta, find omega such that 
 * the refractive index n = k*c/omega matches the cold plasma dispersion.
 * 
 * We use a numerical root-finding approach (Bisection/Brent) on the function:
 * f(omega) = n^2(omega, theta) - (k*c/omega)^2 = 0
 */
real_t WhistlerDispersion::calculateOmega(real_t k, real_t theta) const {
    if (use_simple_dispersion) {
        return calculateOmegaSimple(k, theta);
    }
    
    // Target n^2 as a function of omega
    auto target_n2 = [&](real_t omega) -> real_t {
        real_t kc_omega = (k * c) / omega;
        return kc_omega * kc_omega;
    };
    
    // Residual function: difference between plasma n^2 and target n^2
    // We specifically look for the Whistler branch (Right-hand polarized, usually the larger n^2 root below omega_ce)
    auto residual = [&](real_t omega) -> real_t {
        if (omega <= 0 || omega >= omega_ce) return 1e10; // Whistler exists only for 0 < omega < omega_ce
        
        std::vector<real_t> n2_roots = solveRefractiveIndex(omega, theta);
        if (n2_roots.empty()) return 1e10;
        
        // Select the whistler branch (typically the larger n^2 root, which corresponds to the R-wave/fast mode below omega_ce)
        real_t n2_plasma = *std::max_element(n2_roots.begin(), n2_roots.end());
        
        return n2_plasma - target_n2(omega);
    };
    
    // Search bounds for omega
    // Whistler waves exist between ion cyclotron and electron cyclotron frequencies
    real_t omega_min = 10.0 * omega_ci; 
    real_t omega_max = 0.99 * omega_ce; // Avoid electron cyclotron resonance singularity
    
    // Check if a root exists in the interval
    real_t f_min = residual(omega_min);
    real_t f_max = residual(omega_max);
    
    if (f_min * f_max > 0) {
        // Fallback: If no sign change, try a denser scan to find a bracket
        const int scan_steps = 100;
        real_t d_omega = (omega_max - omega_min) / scan_steps;
        real_t f_prev = f_min;
        real_t w_prev = omega_min;
        bool found_bracket = false;
        
        for (int i = 1; i <= scan_steps; ++i) {
            real_t w_curr = omega_min + i * d_omega;
            real_t f_curr = residual(w_curr);
            if (f_prev * f_curr < 0) {
                omega_min = w_prev;
                omega_max = w_curr;
                f_min = f_prev;
                f_max = f_curr;
                found_bracket = true;
                break;
            }
            f_prev = f_curr;
            w_prev = w_curr;
        }
        
        if (!found_bracket) {
            // std::cerr << "Warning: No whistler root found for k=" << k << ", theta=" << theta*180/M_PI << " deg" << std::endl;
            return -1.0;
        }
    }
    
    // Bisection method to find the root
    const int max_iter = 60;
    const real_t tol = 1e-8 * omega_ce;
    real_t omega_mid = 0.5 * (omega_min + omega_max);
    
    for (int iter = 0; iter < max_iter; ++iter) {
        omega_mid = 0.5 * (omega_min + omega_max);
        real_t f_mid = residual(omega_mid);
        
        if (std::abs(f_mid) < 1e-6 || (omega_max - omega_min) < tol) {
            break;
        }
        
        if (f_min * f_mid < 0) {
            omega_max = omega_mid;
            f_max = f_mid;
        } else {
            omega_min = omega_mid;
            f_min = f_mid;
        }
    }
    
    return omega_mid;
}

/**
 * Legacy simplified dispersion (kept for fallback/testing)
 */
real_t WhistlerDispersion::calculateOmegaSimple(real_t k, real_t theta) const {
    real_t k_parallel = k * std::cos(theta);
    return k * std::abs(k_parallel) * w_factor;
}

// ... (Keep other utility methods like calculateGroupVelocity, printInfo, etc. as they were)