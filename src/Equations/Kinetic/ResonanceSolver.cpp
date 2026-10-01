/**
 * Implementation of the ResonanceSolver class.
 * 
 * Solves the coupled dispersion relation and resonance condition:
 *   ω(k, θ_k) = k cos(θ_k) v_∥ + n Ω_ce / γ
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
    
    // Define resonance residual function using standard resonance condition:
    // f(p) = omega_wave * gamma(p) - k_parallel * c * p * xi - n * omega_ce
    // This avoids division by gamma and is numerically more stable
    // Standard resonance: ω - k_∥·v_∥ - n·Ω_ce/γ = 0
    // Multiply by γ: ω·γ - k_∥·v·γ·ξ - n·Ω_ce = 0
    // Since v·γ = p·c: ω·γ - k_∥·c·p·ξ - n·Ω_ce = 0
    auto resonance_residual = [&](real_t p) -> real_t {
        real_t gamma = std::sqrt(1.0 + p * p);
        return omega_wave * gamma - k_parallel * c * p * xi - static_cast<real_t>(n) * omega_ce;
    };
    
    // DEBUG: Print parameters for first call (disabled for cleaner output)
    /*
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
    */
    
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
            real_t original_residual = omega_wave - k_parallel * v_parallel_root - static_cast<real_t>(n) * omega_ce / gamma_root;
            real_t relative_error = std::abs(original_residual) / (std::abs(omega_wave) + 1e-10);
            
            if (relative_error < 1e-6 && p_root > 0) {
                resonant_p_values.push_back(p_root);
                
                // Debug output disabled
                /*
                if (!debug_printed) {
                    std::cerr << "  Found root: p=" << p_root << ", original_residual=" << original_residual 
                              << ", rel_err=" << relative_error << std::endl;
                }
                */
            }
        }
        
        f_prev = f_curr;
    }
    
    // Debug flag disabled
    // debug_printed = true;
    
    return resonant_p_values;
}
