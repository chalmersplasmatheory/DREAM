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
    // Derive ion frequencies assuming deuterium plasma (Z=1, n_i = n_e)
    // omega_pi = omega_pe * sqrt(m_e / m_D)
    // omega_ci = omega_ce * (m_e / m_D)
    real_t mass_ratio = m_electron / m_ion;  // m_e / m_D
    omega_pi = omega_pe * std::sqrt(mass_ratio);
    omega_ci = omega_ce * mass_ratio;
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
 * Calculate wave polarization (Ex, Ey, Ez) and dielectric derivative (den).
 * Based on cold plasma dielectric tensor.
 *
 * den = [D(ω+δω) - D(ω)] / (δω · ω) where D is the dispersion function.
 */
void WaveParticleCoupling::calculate_den(
    real_t omega, real_t k, real_t theta_k,
    real_t &Ex, real_t &Ey, real_t &Ez, real_t &den
) const {
    real_t cc = c;
    real_t kpar = k * std::cos(theta_k);
    real_t kperp = k * std::sin(theta_k);
    real_t omegace2 = omega_ce * omega_ce;
    real_t omegaci2 = omega_ci * omega_ci;
    real_t omegape2 = omega_pe * omega_pe;
    real_t omegapi2 = omega_pi * omega_pi;
    real_t omega2 = omega * omega;

    // Ex = 1 (normalized)
    Ex = 1.0;

    // Ey from cold plasma polarization
    real_t n2 = k * k * cc * cc / omega2;  // refractive index squared

    Ey = (omega_ce / omega * omegape2 / (omega2 - omegace2)
        - omega_ci / omega * omegapi2 / (omega2 - omegaci2))
       / (1.0 - omegape2 / (omega2 - omegace2)
              - omegapi2 / (omega2 - omegaci2)
              - n2);

    // Ez from parallel component
    Ez = -kpar * cc / omega
       / (1.0 - omegape2 / omega2 - omegapi2 / omega2
              - kperp * kperp * cc * cc / (omega2))
       * kperp * cc / omega;

    // Denominator D(ω) and numerical derivative ∂D/∂ω
    // D(ω) = [(1+Ey²)·S - 2·Ey·T + Ez²·P] · ω²
    // where S = (R+L)/2, T = (R-L)/2 in Stix notation
    // Polarization (Ey, Ez) is FROZEN at omega: only the Stix parameters are
    // differentiated (paper convention, eq. (12)).
    auto D_func = [&](real_t w) -> real_t {
        real_t w2 = w * w;
        real_t S_w = 1.0 - omegape2 / (w2 - omegace2) - omegapi2 / (w2 - omegaci2);
        real_t T_w = omega_ce / w * omegape2 / (w2 - omegace2)
                   - omega_ci / w * omegapi2 / (w2 - omegaci2);
        real_t P_w = 1.0 - omegape2 / w2 - omegapi2 / w2;

        return ((1.0 + Ey * Ey) * S_w - 2.0 * Ey * T_w + Ez * Ez * P_w) * w2;
    };

    real_t D_omega = D_func(omega);
    real_t D_omega_shifted = D_func(omega + delta_omega);
    den = (D_omega_shifted - D_omega) / (delta_omega * omega);
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
 * Calculate the normalized bracket term for quasi-linear diffusion.
 *
 *   bracket = e²/(4π m_e² c²) · |ψ_norm|² / |v_g - v∥ cos θ_k|
 *
 * where |ψ_norm|² is the normalized coupling matrix element with Ex=1:
 *   |ψ_norm|² = |Ex·(n/z)·Jn(z) + i·Ey·Jn'(z) + (p∥/p⊥)·Ez·Jn(z)|²
 *
 * The actual coupling strength is then: bracket × |E_k|²
 * where |E_k|² = amplitude × m_e c² / (ε₀ × den)
 */
real_t WaveParticleCoupling::calculateCouplingStrength(
    real_t p, real_t xi, real_t k, real_t theta_k, int n,
    const WhistlerDispersion &dispersion
) const {
    real_t omega = dispersion.calculateOmega(k, theta_k);
    if (omega < 0) return 0.0;

    // Get wave polarization and dielectric derivative from cold plasma theory
    real_t Ex, Ey, Ez, den;
    calculate_den(omega, k, theta_k, Ex, Ey, Ez, den);

    // Particle kinematics
    real_t gamma = std::sqrt(1.0 + p * p);
    real_t v_parallel = (p / gamma) * xi * c;      // m/s
    real_t p_perp = p * std::sqrt(std::max(0.0, 1.0 - xi * xi));

    // Bessel function argument
    real_t k_perp = k * std::sin(theta_k);
    real_t rho = calculateGyroradius(p, xi);
    real_t z = k_perp * rho;                        // k⊥ρ (dimensionless)

    real_t bessel_n   = besselJ(n, z);
    real_t bessel_np1 = besselJ(n + 1, z);
    real_t bessel_nm1 = besselJ(n - 1, z);

    // Jn'(z) = [J_{n-1}(z) - J_{n+1}(z)] / 2
    real_t bessel_n_prime = (bessel_nm1 - bessel_np1) / 2.0;

    // Coupling term: n·Jn(z)/z
    real_t n_over_z_Jn = 0.0;
    if (std::abs(z) > 1e-10) {
        n_over_z_Jn = static_cast<real_t>(n) * bessel_n / z;
    } else {
        n_over_z_Jn = (n == 0) ? 0.5 : 0.0;
    }

    // Normalized matrix element (Ex = 1):
    // ψ_norm = Ex·(n/z)·Jn + i·Ey·Jn' + (p∥/p⊥)·Ez·Jn
    real_t ppar_over_pperp = 0.0;
    if (p_perp > 1e-30) {
        ppar_over_pperp = (p * xi) / p_perp;
    }

    // |ψ_norm|² = [Re]² + [Im]²
    //   Re = (n/z)·Jn + (p∥/p⊥)·Ez·Jn
    //   Im = Ey·Jn'
    real_t psi_real = n_over_z_Jn + ppar_over_pperp * Ez * bessel_n;
    real_t psi_imag = Ey * bessel_n_prime;
    real_t psi_norm_sq = psi_real * psi_real + psi_imag * psi_imag;

    // Resonance denominator |∂ω/∂k - v∥ cos θ_k|
    real_t v_g = calculateGroupVelocity(k, theta_k, dispersion);
    real_t denom = std::abs(v_g - v_parallel * std::cos(theta_k));
    if (denom < 1e-10) denom = 1e-10;

    // Prefactor: e²/(4π m_e² c²)
    static constexpr real_t prefactor = M_PI * e_charge * e_charge
                                      / (m_electron * m_electron * c * c);

    return prefactor * psi_norm_sq / denom;
}

/**
 * Calculate both the bracket term and dielectric derivative den.
 * Avoids computing omega, polarization, and Bessel functions twice.
 */
void WaveParticleCoupling::calculateBracketAndDen(
    real_t p, real_t xi, real_t k, real_t theta_k, int n,
    const WhistlerDispersion &dispersion,
    real_t &bracket, real_t &den
) const {
    real_t omega = dispersion.calculateOmega(k, theta_k);
    if (omega < 0) { bracket = 0.0; den = 1.0; return; }

    real_t Ex, Ey, Ez;
    calculate_den(omega, k, theta_k, Ex, Ey, Ez, den);

    real_t gamma = std::sqrt(1.0 + p * p);
    real_t v_parallel = (p / gamma) * xi * c;
    real_t p_perp = p * std::sqrt(std::max(0.0, 1.0 - xi * xi));

    real_t k_perp = k * std::sin(theta_k);
    real_t rho = calculateGyroradius(p, xi);
    real_t z = k_perp * rho;

    real_t bessel_n   = besselJ(n, z);
    real_t bessel_np1 = besselJ(n + 1, z);
    real_t bessel_nm1 = besselJ(n - 1, z);
    real_t bessel_n_prime = (bessel_nm1 - bessel_np1) / 2.0;

    real_t n_over_z_Jn = 0.0;
    if (std::abs(z) > 1e-10)
        n_over_z_Jn = static_cast<real_t>(n) * bessel_n / z;
    else
        n_over_z_Jn = (n == 0) ? 0.5 : 0.0;

    real_t ppar_over_pperp = (p_perp > 1e-30) ? (p * xi) / p_perp : 0.0;

    real_t psi_real = n_over_z_Jn + ppar_over_pperp * Ez * bessel_n;
    real_t psi_imag = Ey * bessel_n_prime;
    real_t psi_norm_sq = psi_real * psi_real + psi_imag * psi_imag;

    real_t v_g = calculateGroupVelocity(k, theta_k, dispersion);
    real_t denom = std::abs(v_g - v_parallel * std::cos(theta_k));
    if (denom < 1e-10) denom = 1e-10;

    static constexpr real_t prefactor = M_PI * e_charge * e_charge
                                      / (m_electron * m_electron * c * c);

    bracket = prefactor * psi_norm_sq / denom;
}
