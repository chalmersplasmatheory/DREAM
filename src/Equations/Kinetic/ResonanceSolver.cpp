/**
 * Implementation of the ResonanceSolver class.
 * 
 * Solves the coupled dispersion relation and resonance condition:
 *   ω(k, θ_k) = k cos(θ_k) v_∥ - n Ω_ce / γ
 * 
 * Uses numerical root-finding to locate resonant wavenumbers.
 */

#include "DREAM/Equations/Kinetic/ResonanceSolver.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace DREAM;

/**
 * Constructor
 */
ResonanceSolver::ResonanceSolver(real_t omega_ce_val)
    : omega_ce(omega_ce_val) {
}

/**
 * Calculate parallel velocity (normalized to c)
 */
real_t ResonanceSolver::calculateVparallel(real_t p, real_t xi) const {
    // v/c = p / sqrt(1 + p^2)
    // v_parallel/c = (v/c) * xi
    real_t gamma = std::sqrt(1.0 + p * p);
    return (p / gamma) * xi;
}

/**
 * Calculate Lorentz factor
 */
real_t ResonanceSolver::calculateGamma(real_t p) const {
    return std::sqrt(1.0 + p * p);
}

/**
 * Resonance residual function
 */
real_t ResonanceSolver::resonanceResidual(
    real_t k, real_t vpar_over_c, real_t gamma,
    real_t cos_theta, int n,
    const WhistlerDispersion &dispersion, real_t theta_k
) const {
    // Calculate wave frequency from dispersion relation
    real_t omega_disp = dispersion.calculateOmega(k, theta_k);
    
    // Calculate resonant frequency from particle parameters
    real_t kpar = k * cos_theta;
    real_t vpar = vpar_over_c * c;  // Convert to m/s
    real_t omega_res = kpar * vpar - static_cast<real_t>(n) * omega_ce / gamma;
    
    // Return residual
    return omega_disp - omega_res;
}

/**
 * Find resonant wavenumber(s) using bracketing and bisection
 */
std::vector<real_t> ResonanceSolver::findResonantK(
    real_t p, real_t xi, real_t theta_k, int n,
    const WhistlerDispersion &dispersion,
    real_t k_min, real_t k_max
) const {
    std::vector<real_t> resonant_k_values;
    
    // Calculate particle parameters
    real_t vpar_over_c = calculateVparallel(p, xi);
    real_t gamma = calculateGamma(p);
    real_t cos_theta = std::cos(theta_k);
    
    // Check if search range is valid
    if (k_min >= k_max) {
        std::cerr << "Warning: Invalid k range [" << k_min << ", " << k_max << "]" << std::endl;
        return resonant_k_values;
    }
    
    // Strategy: Sample the residual function on a grid and look for sign changes
    // Then use bisection to find roots in each interval
    
    const int num_samples = 100;  // Number of sample points
    real_t dk = (k_max - k_min) / static_cast<real_t>(num_samples);
    
    real_t f_prev = resonanceResidual(k_min, vpar_over_c, gamma, cos_theta, n, dispersion, theta_k);
    
    for (int i = 1; i <= num_samples; i++) {
        real_t k_curr = k_min + static_cast<real_t>(i) * dk;
        real_t f_curr = resonanceResidual(k_curr, vpar_over_c, gamma, cos_theta, n, dispersion, theta_k);
        
        // Check for sign change (indicates a root in [k_prev, k_curr])
        if (f_prev * f_curr < 0.0) {
            // Found a bracket, use bisection to find root
            real_t k_left = k_curr - dk;
            real_t k_right = k_curr;
            real_t f_left = f_prev;
            real_t f_right = f_curr;
            
            // Bisection method
            const int max_iter = 50;
            real_t k_root = 0.5 * (k_left + k_right);
            
            for (int iter = 0; iter < max_iter; iter++) {
                k_root = 0.5 * (k_left + k_right);
                real_t f_mid = resonanceResidual(k_root, vpar_over_c, gamma, cos_theta, n, dispersion, theta_k);
                
                // Check convergence
                if (std::abs(f_mid) < 1e-6 || (k_right - k_left) / k_root < 1e-8) {
                    break;
                }
                
                // Update bracket
                if (f_left * f_mid < 0.0) {
                    k_right = k_root;
                    f_right = f_mid;
                } else {
                    k_left = k_root;
                    f_left = f_mid;
                }
            }
            
            // Verify the root is physically meaningful
            real_t omega_at_root = dispersion.calculateOmega(k_root, theta_k);
            if (omega_at_root > 0 && k_root > 0) {
                resonant_k_values.push_back(k_root);
            }
        }
        
        f_prev = f_curr;
    }
    
    return resonant_k_values;
}

/**
 * Find resonant wavenumber(s) using pre-calculated omega values.
 * This avoids repeated PDRF calls by using cached dispersion data.
 */
std::vector<real_t> ResonanceSolver::findResonantKWithCachedOmega(
    real_t p, real_t xi, real_t theta_k, int n,
    const real_t* omega_cache,
    const real_t* k_values,
    len_t num_modes
) const {
    std::vector<real_t> resonant_k_values;
    
    // Calculate particle parameters
    real_t vpar_over_c = calculateVparallel(p, xi);
    real_t gamma = calculateGamma(p);
    real_t cos_theta = std::cos(theta_k);
    real_t c_val = 2.99792458e8;  // Speed of light
    
    // Calculate resonant frequency for this particle state
    // We'll check which modes are near resonance
    real_t vpar = vpar_over_c * c_val;
    
    // Loop through all cached modes and find those near resonance
    for (len_t m = 0; m < num_modes; m++) {
        real_t omega_wave = omega_cache[m];
        if (omega_wave < 0) continue;  // Skip invalid modes
        
        real_t k = k_values[m];
        real_t kpar = k * cos_theta;
        
        // Calculate resonance condition residual
        real_t omega_res = kpar * vpar - static_cast<real_t>(n) * omega_ce / gamma;
        real_t residual = omega_wave - omega_res;
        
        // Check if near resonance (within 10% relative error)
        real_t relative_error = std::abs(residual) / (std::abs(omega_wave) + 1e-10);
        if (relative_error < 0.1) {
            resonant_k_values.push_back(k);
        }
    }
    
    return resonant_k_values;
}

/**
 * Check if particle-wave pair is in resonance
 */
bool ResonanceSolver::isInResonance(
    real_t omega, real_t kpar, real_t vpar,
    real_t gamma, int n, real_t tolerance
) const {
    // Calculate resonance condition residual
    real_t residual = omega - kpar * vpar + static_cast<real_t>(n) * omega_ce / gamma;
    
    // Check relative error
    real_t relative_error = std::abs(residual) / (std::abs(omega) + 1e-10);
    
    return relative_error < tolerance;
}

/**
 * Find resonant momentum p for given wave parameters and pitch angle.
 * This uses the REVERSE approach: solve for exact resonant p given (k, theta_k, n, xi).
 * 
 * Uses NUMERICAL root-finding (bisection method) instead of analytical quadratic solution
 * to avoid spurious negative roots introduced by squaring.
 */
std::vector<real_t> ResonanceSolver::findResonantP(
    real_t k, real_t theta_k, int n, real_t xi,
    const WhistlerDispersion &dispersion,
    real_t p_min, real_t p_max
) const {
    std::vector<real_t> resonant_p_values;
    
    // Calculate wave parameters
    real_t omega_wave = dispersion.calculateOmega(k, theta_k);
    real_t cos_theta = std::cos(theta_k);
    real_t k_parallel = k * cos_theta;
    
    // Define resonance residual function using QUADRE-style formulation:
    // f(p) = omega_wave * gamma(p) - k_parallel * c * p * xi + n * omega_ce
    // This avoids division by gamma and is numerically more stable
    // Derived from: omega - k_parallel * v_parallel + n * omega_ce / gamma = 0
    // Multiply by gamma: omega * gamma - k_parallel * (p/gamma) * xi * c * gamma + n * omega_ce = 0
    // Simplify: omega * gamma - k_parallel * c * p * xi + n * omega_ce = 0
    auto resonance_residual = [&](real_t p) -> real_t {
        real_t gamma = std::sqrt(1.0 + p * p);
        return omega_wave * gamma - k_parallel * c * p * xi + static_cast<real_t>(n) * omega_ce;
    };
    
    // DEBUG: Print parameters for first call
    static bool debug_printed = false;
    if (!debug_printed) {
        std::cerr << "DEBUG findResonantP (numerical): k=" << k << ", theta_k=" << theta_k 
                  << ", n=" << n << ", xi=" << xi << std::endl;
        std::cerr << "  omega_wave=" << omega_wave << ", k_parallel=" << k_parallel << std::endl;
        std::cerr << "  Searching in p range [" << p_min << ", " << p_max << "]" << std::endl;
        
        // Sample the residual function to find sign changes
        const int num_samples = 100;
        real_t dp_sample = (p_max - p_min) / num_samples;
        real_t f_prev = resonance_residual(p_min);
        
        std::cerr << "  Sampling residual: f(" << p_min << ")=" << f_prev;
        
        int sign_changes = 0;
        for (int i = 1; i <= num_samples; i++) {
            real_t p_curr = p_min + i * dp_sample;
            real_t f_curr = resonance_residual(p_curr);
            
            if (f_prev * f_curr < 0) {
                sign_changes++;
                if (sign_changes <= 3) {
                    std::cerr << " [sign change near p=" << p_curr << "]";
                }
            }
            f_prev = f_curr;
        }
        std::cerr << " - Total sign changes: " << sign_changes << std::endl;
    }
    
    // Use bisection method to find roots
    // Sample the function on a grid and look for sign changes
    const int num_samples = 200;
    real_t dp = (p_max - p_min) / num_samples;
    
    real_t f_prev = resonance_residual(p_min);
    
    for (int i = 1; i <= num_samples; i++) {
        real_t p_curr = p_min + i * dp;
        real_t f_curr = resonance_residual(p_curr);
        
        // Check for sign change (indicates a root in [p_prev, p_curr])
        if (f_prev * f_curr < 0.0) {
            // Found a bracket, use bisection to find root
            real_t p_left = p_curr - dp;
            real_t p_right = p_curr;
            real_t f_left = f_prev;
            real_t f_right = f_curr;
            
            // Bisection method
            const int max_iter = 50;
            real_t p_root = 0.5 * (p_left + p_right);
            
            for (int iter = 0; iter < max_iter; iter++) {
                p_root = 0.5 * (p_left + p_right);
                real_t f_mid = resonance_residual(p_root);
                
                // Check convergence
                if (std::abs(f_mid) < 1e-6 || (p_right - p_left) / p_root < 1e-8) {
                    break;
                }
                
                // Update bracket
                if (f_left * f_mid < 0.0) {
                    p_right = p_root;
                    f_right = f_mid;
                } else {
                    p_left = p_root;
                    f_left = f_mid;
                }
            }
            
            // Verify the root using original resonance condition
            real_t gamma_root = std::sqrt(1.0 + p_root * p_root);
            real_t v_parallel_root = (p_root / gamma_root) * xi * c;
            real_t original_residual = omega_wave - k_parallel * v_parallel_root + static_cast<real_t>(n) * omega_ce / gamma_root;
            real_t relative_error = std::abs(original_residual) / (std::abs(omega_wave) + 1e-10);
            
            if (relative_error < 1e-6 && p_root > 0) {
                resonant_p_values.push_back(p_root);
                
                if (!debug_printed) {
                    std::cerr << "  Found root: p=" << p_root << ", original_residual=" << original_residual 
                              << ", rel_err=" << relative_error << std::endl;
                }
            }
        }
        
        f_prev = f_curr;
    }
    
    debug_printed = true;
    
    return resonant_p_values;
}
