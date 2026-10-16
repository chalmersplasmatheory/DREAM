/**
 * Implementation of WaveParticleCoupling class.
 * 
 * Calculates the wave-particle coupling strength |Ψ_{n,k}|² based on
 * quasi-linear theory and QUADRE implementation.
 */

#include "DREAM/Equations/Kinetic/WaveParticleCoupling.hpp"
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>  // GSL Bessel functions

using namespace DREAM;

/**
 * Constructor
 */
WaveParticleCoupling::WaveParticleCoupling(
    real_t omega_pe_val, real_t omega_ce_val, real_t n_e_val
) : omega_pe(omega_pe_val), omega_ce(omega_ce_val), n_e(n_e_val) {
}

/**
 * Calculate Bessel function J_n(x) using GSL
 */
real_t WaveParticleCoupling::besselJ(int n, real_t x) const {
    if (x == 0.0) {
        return (n == 0) ? 1.0 : 0.0;
    }
    
    // Use GSL for accurate Bessel function calculation
    if (n >= 0) {
        return gsl_sf_bessel_Jn(n, std::abs(x));
    } else {
        // J_{-n}(x) = (-1)^n J_n(x)
        real_t result = gsl_sf_bessel_Jn(-n, std::abs(x));
        return ((-n) % 2 == 0) ? result : -result;
    }
}

/**
 * Calculate gyroradius ρ = p_⊥/(m_e Ω_ce)
 */
real_t WaveParticleCoupling::calculateGyroradius(real_t p, real_t xi) const {
    // p_⊥ = p · √(1 - ξ²)
    real_t p_perp = p * std::sqrt(std::max(0.0, 1.0 - xi*xi));
    
    // ρ = p_⊥ · m_e · c / (e · B) = p_⊥ · c / Ω_ce
    return p_perp * c / omega_ce;
}

/**
 * Calculate group velocity ∂ω/∂k using numerical differentiation
 */
real_t WaveParticleCoupling::calculateGroupVelocity(
    real_t k, real_t theta_k,
    const WhistlerDispersion &dispersion,
    real_t dk
) const {
    // Central difference: dω/dk ≈ [ω(k+dk) - ω(k-dk)] / (2·dk)
    real_t k_plus = k + dk;
    real_t k_minus = std::max(k - dk, 0.1);  // Avoid k <= 0
    
    real_t omega_plus = dispersion.calculateOmega(k_plus, theta_k);
    real_t omega_minus = dispersion.calculateOmega(k_minus, theta_k);
    
    return (omega_plus - omega_minus) / (2.0 * dk);
}

/**
 * Calculate coupling strength |Ψ_{n,k}|²
 * 
 * Implements the formula from QUADRE inject.py (lines 689-692):
 *   weight = [coupling_term]² / |dω/dk - α| / den · k²·ω_pe²/(n_e·π)
 */
real_t WaveParticleCoupling::calculateCouplingStrength(
    real_t p, real_t xi, real_t k, real_t theta_k, int n,
    const WhistlerDispersion &dispersion
) const {
    // Calculate particle parameters
    real_t gamma = std::sqrt(1.0 + p * p);
    real_t v_parallel_over_c = (p / gamma) * xi;
    real_t alpha = v_parallel_over_c * c * std::cos(theta_k);  // v_∥ cos(θ_k)
    
    // Calculate wave parameters
    real_t omega = dispersion.calculateOmega(k, theta_k);
    real_t k_perp = k * std::sin(theta_k);
    
    // Calculate gyroradius and argument for Bessel functions
    real_t rho = calculateGyroradius(p, xi);
    real_t kperp_rho = k_perp * rho;
    
    // Get wave polarization from dispersion solver
    // For now, use simplified whistler wave polarization
    // In full implementation, this should come from eigenvector calculation
    real_t Ex = 1.0;  // Normalized
    real_t Ey = 0.0;  // Simplified - should be calculated from dispersion
    real_t Ez = 0.0;  // Simplified - should be calculated from dispersion
    
    // TODO: Implement proper polarization calculation
    // This requires solving for eigenvectors of the dispersion matrix
    // For whistler waves, typical polarization is:
    // E_y/E_x ≈ i·ω/Ω_ce (circular polarization)
    // E_z is small for parallel propagation
    
    // Calculate coupling term (from QUADRE line 689)
    real_t bessel_n = besselJ(n, kperp_rho);
    real_t bessel_np1 = besselJ(n + 1, kperp_rho);
    real_t bessel_nm1 = besselJ(n - 1, kperp_rho);
    
    real_t coupling_term = 0.0;
    
    if (std::abs(kperp_rho) > 1e-10) {
        coupling_term += static_cast<real_t>(n) * bessel_n / kperp_rho;
    } else {
        // Limit as kperp_rho → 0
        coupling_term += (n == 0) ? 0.5 : 0.0;
    }
    
    // Add E_z term (small for whistler waves)
    if (std::abs(xi) < 1.0) {  // Avoid division by zero at ξ=±1
        coupling_term += Ez * bessel_n * xi / std::sqrt(1.0 - xi*xi);
    }
    
    // Add E_y term (important for circular polarization)
    coupling_term -= Ey * (bessel_np1 - bessel_nm1) / 2.0;
    
    // Square the coupling term
    real_t weight = coupling_term * coupling_term;
    
    // Calculate denominator |∂ω/∂k - α|
    real_t group_velocity = calculateGroupVelocity(k, theta_k, dispersion);
    real_t denom = std::abs(group_velocity - alpha);
    
    if (denom < 1e-10) {
        // Avoid division by zero
        denom = 1e-10;
    }
    
    weight /= denom;
    
    // Normalize by wave energy density (simplified)
    // Full implementation needs dielectric tensor derivative "den"
    // For now, use simplified normalization
    real_t normalization = k * k * omega_pe * omega_pe / (n_e * M_PI);
    weight *= normalization;
    
    return weight;
}
