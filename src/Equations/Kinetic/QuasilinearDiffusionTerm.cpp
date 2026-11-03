/**
 * Implementation of QuasilinearDiffusionTerm.
 * 
 * Implements quasi-linear diffusion from wave-particle interactions
 * based on QUADRE implementation, converted from FEM to FVM.
 */

#include "DREAM/Equations/Kinetic/QuasilinearDiffusionTerm.hpp"
#include "DREAM/DREAMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include <iostream>
#include <cmath>

using namespace DREAM;

/**
 * Constructor
 */
QuasilinearDiffusionTerm::QuasilinearDiffusionTerm(
    FVM::Grid *grid,
    WaveSpectrum *spectrum,
    WhistlerDispersion *dispersion,
    ResonanceSolver *resonanceSolver,
    WaveParticleCoupling *coupling,
    const std::vector<int> &harmonicModes,
    real_t start_inject_time,
    real_t inject_cycle_duration,
    real_t ramp_time
) : FVM::DiffusionTerm(grid),
    spectrum(spectrum),
    dispersion(dispersion),
    resonanceSolver(resonanceSolver),
    coupling(coupling),
    harmonicModes(harmonicModes),
    omega_cache(nullptr),
    dispersion_cached(false),
    D_pp_cache(nullptr),
    D_pxi_pface_cache(nullptr),
    D_pxi_xface_cache(nullptr),
    D_xixi_cache(nullptr),
    nr_cached(0),
    np1_cached(0),
    np2_cached(0),
    first_rebuild_done(false),
    operator_cached(false),
    last_amplitudes(nullptr),
    num_modes_tracked(0),
    start_inject_time(start_inject_time),
    inject_cycle_duration(inject_cycle_duration),
    ramp_time(ramp_time),
    base_amplitude(spectrum ? spectrum->getAmplitude(0) : 1.0) {
    
    SetName("QuasilinearDiffusionTerm");
    
    // Pre-calculate and cache dispersion relation for all modes
    precalculateDispersion();
}

/**
 * Destructor - free cached arrays
 */
QuasilinearDiffusionTerm::~QuasilinearDiffusionTerm() {
    if (omega_cache) delete[] omega_cache;
    if (D_pp_cache) delete[] D_pp_cache;
    if (D_pxi_pface_cache) delete[] D_pxi_pface_cache;
    if (D_pxi_xface_cache) delete[] D_pxi_xface_cache;
    if (D_xixi_cache) delete[] D_xixi_cache;
    if (last_amplitudes) delete[] last_amplitudes;
}

/**
 * Calculate effective wave amplitude at time t based on periodic injection schedule.
 * 
 * Features:
 * - After start_inject_time / at each ON cycle start: linear ramp-up over ramp_time
 * - OFF transitions: immediate step to 0
 * - 50% duty cycle (first half ON, second half OFF)
 */
real_t QuasilinearDiffusionTerm::calculateEffectiveAmplitude(real_t t) const {
    // Case 1: Injection hasn't started yet
    if (start_inject_time >= 0 && t < start_inject_time) {
        return 0.0;
    }
    
    // Time since injection started (or 0 if immediate)
    real_t t_local = (start_inject_time >= 0) ? (t - start_inject_time) : t;
    
    // Case 2: Continuous injection (no periodicity)
    if (inject_cycle_duration <= 0) {
        if (ramp_time > 0 && t_local < ramp_time)
            return base_amplitude * (t_local / ramp_time);
        return base_amplitude;
    }
    
    // Case 3: Periodic injection with 50% duty cycle
    real_t time_in_cycle = std::fmod(t_local, inject_cycle_duration);
    
    // First half of cycle: injection ON (with ramp-up at start)
    if (time_in_cycle < inject_cycle_duration / 2.0) {
        if (ramp_time > 0 && time_in_cycle < ramp_time)
            return base_amplitude * (time_in_cycle / ramp_time);
        return base_amplitude;
    }
    // Second half of cycle: injection OFF (immediate drop)
    else {
        return 0.0;
    }
}

/**
 * Get cached D_pp coefficient
 */
real_t QuasilinearDiffusionTerm::getD_pp(len_t ir, len_t i, len_t j) const {
    if (!D_pp_cache || ir >= nr_cached) return 0.0;
    len_t idx = ir * (np1_cached + 1) * np2_cached + i * np2_cached + j;
     return D_pp_cache[idx];
}

/**
 * Get cached D_pξ coefficient (ξ-face layout, for backward compatibility)
 */
real_t QuasilinearDiffusionTerm::getD_pxi(len_t ir, len_t i, len_t j) const {
    if (!D_pxi_xface_cache || ir >= nr_cached) return 0.0;
    len_t idx = ir * np1_cached * (np2_cached + 1) + i * (np2_cached + 1) + j;
    return D_pxi_xface_cache[idx];
}

/**
 * Get cached D_ξξ coefficient
 */
real_t QuasilinearDiffusionTerm::getD_xixi(len_t ir, len_t i, len_t j) const {
    if (!D_xixi_cache || ir >= nr_cached) return 0.0;
    len_t idx = ir * np1_cached * (np2_cached + 1) + i * (np2_cached + 1) + j;
    return D_xixi_cache[idx];
}

/**
 * Pre-calculate and cache dispersion relation for all wave modes.
 * This is done once at initialization to avoid repeated PDRF calls.
 */
void QuasilinearDiffusionTerm::precalculateDispersion() {
    const len_t numModes = spectrum->getNumModes();
    
    if (omega_cache) delete[] omega_cache;
    omega_cache = new real_t[numModes];
    
    std::cerr << "Pre-calculating dispersion relation for " << numModes << " wave modes..." << std::endl;
    
    for (len_t m = 0; m < numModes; m++) {
        real_t k = spectrum->getK(m);
        real_t theta_k = spectrum->getKtheta(m);
        
        omega_cache[m] = dispersion->calculateOmega(k, theta_k);
        
        // Display progress every 10% or at completion
        if ((m + 1) % (numModes / 10) == 0 || m == numModes - 1) {
            int percent = static_cast<int>((m + 1) * 100.0 / numModes);
            std::cerr << "  Dispersion: " << percent << "% complete (" << (m + 1) << "/" << numModes << " modes)" << std::endl;
        }
    }
    
    dispersion_cached = true;
    std::cerr << "✓ Dispersion relation cached successfully!" << std::endl;
}

/**
 * Rebuild diffusion term coefficients
 * 
 * IMPORTANT: The quasi-linear diffusion operator depends only on:
 *   1. Wave parameters (k, θ_k, ω) - fixed
 *   2. Plasma parameters (n_e, B0) - fixed  
 *   3. Momentum grid geometry (shape functions) - fixed
 * 
 * It does NOT depend on the distribution function f(p, ξ).
 * Therefore, we calculate it once and cache it for reuse.
 */
void QuasilinearDiffusionTerm::Rebuild(const real_t t, const real_t /*dt*/, FVM::UnknownQuantityHandler* /*uqh*/) {
    // ====================================================================
    // Update wave amplitude based on periodic injection schedule
    // ====================================================================
    real_t effective_amplitude = calculateEffectiveAmplitude(t);
    
    // Update spectrum amplitudes
    if (spectrum) {
        len_t numModes = spectrum->getNumModes();
        bool amplitude_changed = false;
        real_t prev_amp = spectrum->getAmplitude(0);
        
        for (len_t m = 0; m < numModes; m++) {
            real_t old_amp = spectrum->getAmplitude(m);
            if (std::abs(old_amp - effective_amplitude) > 1e-30) {
                spectrum->setAmplitude(m, effective_amplitude);
                amplitude_changed = true;
            }
        }
        
        if (amplitude_changed) {
            // Report only the first update and wave on/off transitions
            // (quiet during ramp steps to avoid one line per time step)
            bool crossed = (prev_amp > 1e-30) != (effective_amplitude > 1e-30);
            if (!first_rebuild_done || crossed)
                std::cerr << "Wave amplitude updated at t=" << t 
                          << " s: A = " << effective_amplitude << std::endl;
            // NOTE: no operator re-assembly here. The diffusion operator is
            // LINEAR in the wave amplitude (E_k^2 ~ amplitude), so the cached
            // operator is assembled once at a reference amplitude and scaled
            // by amp_scale when filling the FVM matrices below.
        }
    }
    
    const len_t nr = grid->GetNr();
    
    // Allocate or reallocate cache if needed (first time or grid size changed)
    if (!operator_cached || nr != nr_cached) {
        if (D_pp_cache) delete[] D_pp_cache;
        if (D_pxi_pface_cache) delete[] D_pxi_pface_cache;
        if (D_pxi_xface_cache) delete[] D_pxi_xface_cache;
        if (D_xixi_cache) delete[] D_xixi_cache;
        
        auto *mg = grid->GetMomentumGrid(0);
        np1_cached = mg->GetNp1();
        np2_cached = mg->GetNp2();
        
        D_pp_cache        = new real_t[nr * (np1_cached + 1) * np2_cached];
        D_pxi_pface_cache = new real_t[nr * (np1_cached + 1) * np2_cached];
        D_pxi_xface_cache = new real_t[nr * np1_cached * (np2_cached + 1)];
        D_xixi_cache      = new real_t[nr * np1_cached * (np2_cached + 1)];
        
        nr_cached = nr;
    }
    
    // Calculate or load diffusion coefficients for each radial point (only if not cached).
    // The operator is assembled ONCE at a fixed reference amplitude (linear in A),
    // then scaled per time step by amp_scale = A(t)/A_ref.
    if (!operator_cached) {
        real_t A_ref = (base_amplitude > 0) ? base_amplitude : 1.0;
        // Temporarily set the spectrum to the reference amplitude for assembly
        if (spectrum) {
            len_t numModes = spectrum->getNumModes();
            for (len_t m = 0; m < numModes; m++)
                spectrum->setAmplitude(m, A_ref);
        }

        // Use REVERSE resonance solving method for grid-independent results
        for (len_t ir = 0; ir < nr; ir++) {
            calculateDiffusionCoefficientsReverse(ir);
        }

        // Restore the physical amplitude on the spectrum
        if (spectrum) {
            len_t numModes = spectrum->getNumModes();
            for (len_t m = 0; m < numModes; m++)
                spectrum->setAmplitude(m, effective_amplitude);
        }
        operator_ref_amplitude = A_ref;

        if (!first_rebuild_done)
            std::cerr << "✓ Quasi-linear diffusion operator computed on-the-fly (REVERSE method)" << std::endl;
    }
    
    // Mark operator as cached - no need to recalculate in future time steps
    operator_cached = true;

    // Scale factor from the cached (reference-amplitude) operator to the
    // current amplitude: D is exactly linear in the wave amplitude.
    real_t amp_scale = (operator_ref_amplitude > 0)
                       ? effective_amplitude / operator_ref_amplitude
                       : 0.0;
    
    // Fill FVM diffusion matrices
    // Similar to SynchrotronTerm pattern but using calculated D coefficients
    
    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        
        // Fill diffusion coefficient matrices with NaN check
        // D11 = D_pp (p-p diffusion at f1 grid points)
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1 + 1; i++) {
                real_t D_pp_val = amp_scale * getD_pp(ir, i, j);
                
                // Final safety check before filling FVM matrix
                if (std::isnan(D_pp_val) || std::isinf(D_pp_val)) {
                    std::cerr << "FATAL: Attempting to fill NaN/Inf into D11 at ir=" << ir 
                             << ", i=" << i << ", j=" << j << std::endl;
                    throw DREAM::DREAMException("NaN/Inf detected in D_pp during FVM matrix assembly");
                }
                
                D11(ir, i, j) += D_pp_val;
            }
        }
        
        // D22 = D_ξξ (ξ-ξ diffusion at f2 grid points)
        for (len_t j = 0; j < np2 + 1; j++) {
            for (len_t i = 0; i < np1; i++) {
                real_t D_xixi_val = amp_scale * getD_xixi(ir, i, j);
                
                // Final safety check
                if (std::isnan(D_xixi_val) || std::isinf(D_xixi_val)) {
                    std::cerr << "FATAL: Attempting to fill NaN/Inf into D22 at ir=" << ir 
                             << ", i=" << i << ", j=" << j << std::endl;
                    throw DREAM::DREAMException("NaN/Inf detected in D_xixi during FVM matrix assembly");
                }
                
                D22(ir, i, j) += D_xixi_val;
            }
        }
        
        // D12 = D_pξ at p-face (same grid as D11: (np1+1) × np2)
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1 + 1; i++) {
                real_t val = amp_scale * D_pxi_pface_cache[ir * (np1 + 1) * np2 + i * np2 + j];
                if (std::isnan(val) || std::isinf(val)) {
                    std::cerr << "FATAL: NaN/Inf in D12 (p-face) at ir=" << ir
                             << ", i=" << i << ", j=" << j << std::endl;
                    throw DREAM::DREAMException("NaN/Inf in D_pxi p-face during FVM assembly");
                }
                D12(ir, i, j) += val;
            }
        }
        
        // D21 = D_pξ at ξ-face (same grid as D22: np1 × (np2+1))
        for (len_t j = 0; j < np2 + 1; j++) {
            for (len_t i = 0; i < np1; i++) {
                real_t val = amp_scale * D_pxi_xface_cache[ir * np1 * (np2 + 1) + i * (np2 + 1) + j];
                if (std::isnan(val) || std::isinf(val)) {
                    std::cerr << "FATAL: NaN/Inf in D21 (ξ-face) at ir=" << ir
                             << ", i=" << i << ", j=" << j << std::endl;
                    throw DREAM::DREAMException("NaN/Inf in D_pxi ξ-face during FVM assembly");
                }
                D21(ir, i, j) += val;
            }
        }
        
        // ===== VALIDATION: Print diffusion coefficient stats =====
        // NOTE: This must be done AFTER filling the FVM matrices.
        // Print only on wave on/off transitions (or the very first build).
        // Caches hold the reference-amplitude operator; stats are scaled to
        // the current amplitude.
        bool wave_on_now = (spectrum ? spectrum->getAmplitude(0) : 0.0) > 1e-30;
        bool state_changed = (validation_state == -1) || ((validation_state == 1) != wave_on_now);
        if (ir == 0 && state_changed) {
            validation_state = wave_on_now ? 1 : 0;
            auto *mg = grid->GetMomentumGrid(ir);
            const len_t np1 = mg->GetNp1();
            const len_t np2 = mg->GetNp2();
            
            real_t max_D11 = 0.0, max_D22 = 0.0, max_D12 = 0.0;
            real_t sum_D11 = 0.0, sum_D22 = 0.0, sum_D12 = 0.0;
            len_t count_nonzero_D11 = 0, count_nonzero_D22 = 0, count_nonzero_D12 = 0;
            
            // D11: p-face grid (np1+1) × np2 — read directly from cache
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1 + 1; i++) {
                    real_t val = std::abs(D_pp_cache[ir*(np1+1)*np2 + i*np2 + j]);
                    if (val > max_D11) max_D11 = val;
                    sum_D11 += val;
                    if (val > 1e-30) count_nonzero_D11++;
                }
            }
            
            // D22: xi-face grid np1 × (np2+1) — read directly from cache
            for (len_t j = 0; j < np2 + 1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    real_t val = std::abs(D_xixi_cache[ir*np1*(np2+1) + i*(np2+1) + j]);
                    if (val > max_D22) max_D22 = val;
                    sum_D22 += val;
                    if (val > 1e-30) count_nonzero_D22++;
                }
            }
            
            // D12 p-face: (np1+1) × np2 — independent loop, includes boundary i=np1
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1 + 1; i++) {
                    real_t val = std::abs(D_pxi_pface_cache[ir*(np1+1)*np2 + i*np2 + j]);
                    if (val > max_D12) max_D12 = val;
                    sum_D12 += val;
                    if (val > 1e-30) count_nonzero_D12++;
                }
            }
            
            std::cerr << "\n=== QUASILINEAR DIFFUSION VALIDATION (r=1/" << grid->GetNr() << ") ===" << std::endl;
            std::cerr << "D11 (p-p): max=" << max_D11 * amp_scale << ", sum=" << sum_D11 * amp_scale
                      << ", nonzero=" << count_nonzero_D11 << "/" << ((np1+1)*np2) << std::endl;
            std::cerr << "D22 (ξ-ξ): max=" << max_D22 * amp_scale << ", sum=" << sum_D22 * amp_scale
                      << ", nonzero=" << count_nonzero_D22 << "/" << (np1*(np2+1)) << std::endl;
            std::cerr << "D12 (p-ξ): max=" << max_D12 * amp_scale << ", sum=" << sum_D12 * amp_scale
                      << ", nonzero=" << count_nonzero_D12 << "/" << ((np1+1)*np2) << std::endl;
            
            if (amp_scale * max_D11 < 1e-30 && amp_scale * max_D22 < 1e-30 && amp_scale * max_D12 < 1e-30) {
                std::cerr << "WARNING: All diffusion coefficients are essentially zero!" << std::endl;
                std::cerr << "  Possible causes:" << std::endl;
                std::cerr << "  1. Wave amplitude is too small" << std::endl;
                std::cerr << "  2. Resonance condition not satisfied in the momentum grid" << std::endl;
                std::cerr << "  3. Incorrect plasma parameters (B0, n_e, etc.)" << std::endl;
            } else {
                std::cerr << "✓ Quasilinear diffusion operator is NON-ZERO and active!" << std::endl;
                
                // Additional check: Verify FVM matrix values are finite
                bool fvm_matrix_ok = true;
                for (len_t j = 0; j < np2 && fvm_matrix_ok; j++) {
                    for (len_t i = 0; i < np1 + 1 && fvm_matrix_ok; i++) {
                        if (std::isnan(D11(ir, i, j)) || std::isinf(D11(ir, i, j))) {
                            std::cerr << "ERROR: D11 FVM matrix has NaN/Inf at i=" << i << ", j=" << j << std::endl;
                            fvm_matrix_ok = false;
                        }
                    }
                }
                for (len_t j = 0; j < np2 + 1 && fvm_matrix_ok; j++) {
                    for (len_t i = 0; i < np1 && fvm_matrix_ok; i++) {
                        if (std::isnan(D22(ir, i, j)) || std::isinf(D22(ir, i, j))) {
                            std::cerr << "ERROR: D22 FVM matrix has NaN/Inf at i=" << i << ", j=" << j << std::endl;
                            fvm_matrix_ok = false;
                        }
                        if (j < np2) {
                            if (std::isnan(D12(ir, i, j)) || std::isinf(D12(ir, i, j))) {
                                std::cerr << "ERROR: D12 FVM matrix has NaN/Inf at i=" << i << ", j=" << j << std::endl;
                                fvm_matrix_ok = false;
                            }
                        }
                    }
                }
                if (fvm_matrix_ok) {
                    std::cerr << "✓ FVM diffusion matrices validated (all finite)" << std::endl;
                } else {
                    std::cerr << "✗ FVM diffusion matrices contain invalid values!" << std::endl;
                }
            }
            std::cerr << "==========================================================\n" << std::endl;
        }
        first_rebuild_done = true;
    }
}

/**
 * Calculate diffusion coefficients using REVERSE resonance solving approach.
 * 
 * Instead of iterating over grid points and checking if they resonate,
 * this method:
 * 1. For each wave mode (k, theta_k, n), solves for exact resonant p given xi
 * 2. Finds which grid cells contain the resonant point
 * 3. Distributes contribution to surrounding grid points using triangular weights
 * 
 * This ensures grid-independent results and captures resonance even when
 * the exact resonant point falls between grid nodes.
 */
void QuasilinearDiffusionTerm::calculateDiffusionCoefficientsReverse(len_t ir) {
    auto *mg = grid->GetMomentumGrid(ir);
    const len_t np1 = mg->GetNp1();
    const len_t np2 = mg->GetNp2();
    
    // Safety check: ensure grid dimensions match cached dimensions
    if (np1 != np1_cached || np2 != np2_cached) {
        std::cerr << "WARNING: Radial point " << ir << " has different grid dimensions (" 
                  << np1 << "x" << np2 << ") than cached (" 
                  << np1_cached << "x" << np2_cached << "). Skipping." << std::endl;
        return;
    }
    
    const len_t numModes = spectrum->getNumModes();
    const len_t numHarmonics = harmonicModes.size();
    
    // Display progress (first assembly only; recalculations are silent)
    if (!first_rebuild_done)
        std::cerr << "  Assembling diffusion operator (REVERSE method) at r=" << ir+1 << "/" << grid->GetNr() 
                  << " (" << numModes << " modes × " << numHarmonics << " harmonics)..." << std::endl;
    
    // Initialize caches to zero
    for (len_t i = 0; i < (np1 + 1) * np2; i++) {
        D_pp_cache[ir * (np1 + 1) * np2 + i] = 0.0;
        D_pxi_pface_cache[ir * (np1 + 1) * np2 + i] = 0.0;
    }
    
    for (len_t i = 0; i < np1 * (np2 + 1); i++) {
        D_pxi_xface_cache[ir * np1 * (np2 + 1) + i] = 0.0;
        D_xixi_cache[ir * np1 * (np2 + 1) + i] = 0.0;
    }
    
    // Resonance width parameter (fraction of grid spacing)
    // This controls how wide the resonance peak is spread
    const real_t resonance_width_factor = 0.5;  // 50% of grid spacing
    
    // Progress tracking
    len_t mode_count = 0;
    int total_resonance_points = 0;
    
    // Quadrature weight for ∫ dk d(cos θ_k)
    len_t num_k = spectrum->getNumK();
    len_t num_ktheta = spectrum->getNumKtheta();
    real_t dk = num_k > 1 ? (spectrum->getK(1) - spectrum->getK(0)) : 1.0;
    real_t cos_theta_min = std::cos(spectrum->getKtheta(0));
    real_t cos_theta_max = std::cos(spectrum->getKtheta(num_k * (num_ktheta - 1)));
    real_t dcos_theta = num_ktheta > 1 ? std::abs(cos_theta_max - cos_theta_min) / num_ktheta : 1.0;
    real_t cell_area = 2.0 * M_PI * dk * dcos_theta;   // int d^3k = 2*pi*k^2*dk*dcos(theta) per mode
    
    // Loop over all wave modes
    for (len_t m = 0; m < numModes; m++) {
        real_t k = spectrum->getK(m);
        real_t theta_k = spectrum->getKtheta(m);
        real_t amplitude = spectrum->getAmplitude(m);
        
        if (amplitude < 1e-30) continue;  // Skip negligible modes
        
        mode_count++;
        
        // Display progress every 25% (first assembly only)
        if (!first_rebuild_done && (mode_count % (numModes / 4) == 0 || m == numModes - 1)) {
            int percent = static_cast<int>(mode_count * 100.0 / numModes);
            std::cerr << "    Mode progress: " << percent << "% (" << mode_count << "/" << numModes << ")" << std::endl;
        }
        
        // Use cached omega
        real_t omega = omega_cache[m];
        if (omega < 0) continue;
        
        real_t k_parallel = k * std::cos(theta_k);
        
        // Loop over harmonic numbers
        for (int n : harmonicModes) {
            // Sample pitch angle grid points
            for (len_t j = 0; j < np2; j++) {
                real_t xi = mg->GetP2(j);  // xi at cell center
                real_t xi_abs = std::abs(xi);
                
                // Avoid singularities at ξ = ±1
                if (xi_abs > 0.999) continue;
                
                // ⭐ REVERSE APPROACH: Solve for exact resonant p
                std::vector<real_t> resonant_ps = resonanceSolver->findResonantP(
                    k, theta_k, n, xi, *dispersion, 0.01, 100.0
                );
                
                // DEBUG: Resonance finding (disabled for cleaner output)
                // static std::map<int, int> debug_count_per_harmonic;
                // if (debug_count_per_harmonic.find(n) == debug_count_per_harmonic.end()) {
                //     debug_count_per_harmonic[n] = 0;
                // }
                // 
                // // Only print for xi close to ±1 (within 0.1)
                // bool xi_near_boundary = (std::abs(std::abs(xi) - 1.0) < 0.1);
                // 
                // if (xi_near_boundary && debug_count_per_harmonic[n] < 3 && ir == 0 && m == 0) {
                //     std::cerr << "DEBUG findResonantP: mode=" << m << ", k=" << k << ", theta_k=" << theta_k 
                //               << ", n=" << n << ", xi=" << xi 
                //               << ", found " << resonant_ps.size() << " solutions" << std::endl;
                //     for (auto p_res : resonant_ps) {
                //         std::cerr << "  -> p_res=" << p_res << std::endl;
                //     }
                //     debug_count_per_harmonic[n]++;
                // }
                
                if (resonant_ps.empty()) continue;  // No resonance for this (xi, n)
                
                // Process each resonant momentum
                for (real_t p_res : resonant_ps) {
                    total_resonance_points++;
                    
                    // Find which grid cell contains p_res
                    // p grid is from 0 to np1 (f1 flux faces)
                    len_t i_low = 0;
                    len_t i_high = np1;
                    
                    // Binary search to find the cell
                    while (i_low < i_high - 1) {
                        len_t i_mid = (i_low + i_high) / 2;
                        real_t p_mid = mg->GetP1_f(i_mid);
                        if (p_res < p_mid) {
                            i_high = i_mid;
                        } else {
                            i_low = i_mid;
                        }
                    }
                    
                    // Now p_res is between p[i_low] and p[i_high]
                    real_t p_low = mg->GetP1_f(i_low);
                    real_t p_high = mg->GetP1_f(i_high);
                    real_t dp = p_high - p_low;
                    
                    // Calculate distance from resonant point to grid nodes
                    real_t dist_to_low = std::abs(p_res - p_low);
                    real_t dist_to_high = std::abs(p_res - p_high);
                    
                    // Define resonance width (triangle base)
                    real_t resonance_width = resonance_width_factor * dp;
                    
                    // Calculate triangular weights
                    // Weight is maximum at p_res and decreases linearly to zero at p_res ± width
                    real_t weight_low = 0.0;
                    real_t weight_high = 0.0;
                    
                    if (dist_to_low < resonance_width) {
                        weight_low = 1.0 - dist_to_low / resonance_width;
                    }
                    if (dist_to_high < resonance_width) {
                        weight_high = 1.0 - dist_to_high / resonance_width;
                    }
                    
                    // Normalize weights so they sum to 1
                    real_t weight_sum = weight_low + weight_high;
                    if (weight_sum < 1e-10) continue;  // Too far from resonance
                    
                    weight_low /= weight_sum;
                    weight_high /= weight_sum;
                    
                    // Calculate bracket term and dielectric derivative den
                    real_t bracket = 0.0, den = 1.0;
                    coupling->calculateBracketAndDen(
                        p_res, xi, k, theta_k, n, *dispersion,
                        bracket, den
                    );
                    
                    // Wave electric field energy (paper convention):
                    // |E_kx|^2 = 4 * E_k * m_e c^2 / (eps0 * D), amplitude = E_k,
                    // and D = den (dielectric derivative with frozen polarization)
                    static constexpr real_t m_e = 9.1094e-31;
                    static constexpr real_t c = 2.99792458e8;
                    static constexpr real_t eps0 = 8.854187817e-12;
                    real_t E_k_sq = 4.0 * amplitude * m_e * c * c / (eps0 * den);
                    
                    // Full contribution: bracket × |E_k|² × k² × Δcosθ × weight
                    real_t contrib = bracket * E_k_sq * k * k * cell_area * weight_sum;
                    
                    // Geometric factor at resonant point
                    real_t gamma = std::sqrt(1.0 + p_res * p_res);
                    real_t v_over_c = p_res / gamma;
                    real_t kpar_v_over_omega = k_parallel * v_over_c * c / omega;
                    real_t geom = -xi + kpar_v_over_omega;
                    real_t sqrt_1mxi2 = std::sqrt(1.0 - xi * xi);
                    real_t contrib_pxi   = sqrt_1mxi2 * geom * contrib;
                    real_t contrib_xixi  = geom * geom * contrib;
                    
                    // ---- D_pp: p-face (i_f, j), f1 grid ----
                    len_t idx_pp_low  = ir * (np1 + 1) * np2 + i_low  * np2 + j;
                    len_t idx_pp_high = ir * (np1 + 1) * np2 + i_high * np2 + j;
                    D_pp_cache[idx_pp_low]  += (1.0 - xi * xi) * contrib * weight_low;
                    D_pp_cache[idx_pp_high] += (1.0 - xi * xi) * contrib * weight_high;
                    
                    // ---- D_pξ: p-face (i_f, j), same face grid as D_pp ----
                    // Stencil reads D12 at interior p-faces (i_f < np1, 0 < j < np2-1)
                    len_t idx_pxi_p_low  = ir * (np1 + 1) * np2 + i_low  * np2 + j;
                    len_t idx_pxi_p_high = ir * (np1 + 1) * np2 + i_high * np2 + j;
                    D_pxi_pface_cache[idx_pxi_p_low]  += contrib_pxi * weight_low;
                    D_pxi_pface_cache[idx_pxi_p_high] += contrib_pxi * weight_high;
                    
                    // ---- D_ξξ and D_pξ: ξ-face (i_cell, j_f), f2 grid ----
                    // Resonance at cell center j sits between ξ-faces j and j+1.
                    // p_res is between p-faces i_low and i_high → cell center index = i_low
                    len_t i_cell = i_low;
                    len_t idx_xface_j   = ir * np1 * (np2 + 1) + i_cell * (np2 + 1) + j;
                    len_t idx_xface_jp1 = ir * np1 * (np2 + 1) + i_cell * (np2 + 1) + (j + 1);
                    D_xixi_cache[idx_xface_j]       += 0.5 * contrib_xixi;
                    D_xixi_cache[idx_xface_jp1]     += 0.5 * contrib_xixi;
                    D_pxi_xface_cache[idx_xface_j]  += 0.5 * contrib_pxi;
                    D_pxi_xface_cache[idx_xface_jp1] += 0.5 * contrib_pxi;
                }
            }
        }
    }
    
    if (!first_rebuild_done)
        std::cerr << "  Total resonance points found: " << total_resonance_points << std::endl;
    
    // Save current amplitudes for change detection (same as in calculateDiffusionCoefficients)
    if (!last_amplitudes || num_modes_tracked != numModes) {
        if (last_amplitudes) delete[] last_amplitudes;
        last_amplitudes = new real_t[numModes];
        num_modes_tracked = numModes;
    }
    for (len_t m = 0; m < numModes; m++) {
        last_amplitudes[m] = spectrum->getAmplitude(m);
    }
}
