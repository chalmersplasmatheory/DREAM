/**
 * Implementation of QuasilinearDiffusionTerm.
 * 
 * Implements quasi-linear diffusion from wave-particle interactions
 * based on QUADRE implementation, converted from FEM to FVM.
 */

#include "DREAM/Equations/Kinetic/QuasilinearDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/QLMatrixLoader.hpp"
#include "DREAM/DREAMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include <iostream>
#include <cmath>
#include <map>

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
    const std::vector<int> &harmonicModes
) : FVM::DiffusionTerm(grid),
    spectrum(spectrum),
    dispersion(dispersion),
    resonanceSolver(resonanceSolver),
    coupling(coupling),
    harmonicModes(harmonicModes),
    matrix_loader(nullptr),
    use_precomputed_matrix(false),
    omega_cache(nullptr),
    dispersion_cached(false),
    D_pp_cache(nullptr),
    D_pxi_cache(nullptr),
    D_xixi_cache(nullptr),
    nr_cached(0),
    np1_cached(0),
    np2_cached(0),
    first_rebuild_done(false),
    operator_cached(false),
    current_amplitude(1.0),
    last_amplitudes(nullptr),
    num_modes_tracked(0) {
    
    SetName("QuasilinearDiffusionTerm");
    
    // Pre-calculate and cache dispersion relation for all modes
    precalculateDispersion();
}

/**
 * Constructor using pre-computed matrix from HDF5
 */
QuasilinearDiffusionTerm::QuasilinearDiffusionTerm(
    FVM::Grid *grid,
    const std::string &hdf5_file,
    real_t initial_amplitude
) : FVM::DiffusionTerm(grid),
    spectrum(nullptr),
    dispersion(nullptr),
    resonanceSolver(nullptr),
    coupling(nullptr),
    matrix_loader(nullptr),
    use_precomputed_matrix(true),
    omega_cache(nullptr),
    dispersion_cached(false),
    D_pp_cache(nullptr),
    D_pxi_cache(nullptr),
    D_xixi_cache(nullptr),
    nr_cached(0),
    np1_cached(0),
    np2_cached(0),
    first_rebuild_done(false),
    operator_cached(false),
    current_amplitude(initial_amplitude),
    last_amplitudes(nullptr),
    num_modes_tracked(0) {
    
    SetName("QuasilinearDiffusionTerm");
    
    // Load pre-computed matrix from HDF5
    matrix_loader = new Equations::Kinetic::QLMatrixLoader(hdf5_file);
    
    // Verify grid compatibility
    auto *mg = grid->GetMomentumGrid(0);
    if (mg) {
        len_t np1 = mg->GetNp1();
        len_t np2 = mg->GetNp2();
        
        if (matrix_loader->getNumP() != np1 || matrix_loader->getNumXi() != np2) {
            std::cerr << "Warning: Momentum grid mismatch!" << std::endl;
            std::cerr << "  HDF5 file: " << matrix_loader->getNumP() << " x " << matrix_loader->getNumXi() << std::endl;
            std::cerr << "  DREAM grid: " << np1 << " x " << np2 << std::endl;
            std::cerr << "  Diffusion coefficients may be incorrect." << std::endl;
        }
    }
    
    std::cerr << "✓ Quasi-linear diffusion term initialized with pre-computed matrix" << std::endl;
    std::cerr << "  Initial amplitude: " << initial_amplitude << std::endl;
}

/**
 * Destructor - free cached arrays
 */
QuasilinearDiffusionTerm::~QuasilinearDiffusionTerm() {
    if (omega_cache) delete[] omega_cache;
    if (D_pp_cache) delete[] D_pp_cache;
    if (D_pxi_cache) delete[] D_pxi_cache;
    if (D_xixi_cache) delete[] D_xixi_cache;
    if (last_amplitudes) delete[] last_amplitudes;
    if (matrix_loader) delete matrix_loader;
}

/**
 * Set current wave amplitude (for pre-computed matrix mode)
 */
void QuasilinearDiffusionTerm::setCurrentAmplitude(real_t A_t) {
    current_amplitude = A_t;
    
    if (use_precomputed_matrix && matrix_loader) {
        matrix_loader->setCurrentAmplitude(A_t);
    }
    
    // Mark operator as not cached so it will be rebuilt with new amplitude
    operator_cached = false;
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
 * Get cached D_pξ coefficient
 */
real_t QuasilinearDiffusionTerm::getD_pxi(len_t ir, len_t i, len_t j) const {
    if (!D_pxi_cache || ir >= nr_cached) return 0.0;
    len_t idx = ir * np1_cached * (np2_cached + 1) + i * (np2_cached + 1) + j;
    return D_pxi_cache[idx];
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
 * Quick resonance check: determine if a grid point can possibly resonate
 * with any wave mode. This is a fast pre-filter before calling the full
 * resonance solver.
 * 
 * Based on QUADRE's approach:
 *   omega_res = k*v_parallel*cos(theta_k) + n*omega_ce/gamma
 *   
 * For a given (p, xi), calculate the range of possible resonant frequencies
 * as k varies over [k_min, k_max]. If this range doesn't overlap with
 * the wave spectrum frequency range, skip this point.
 */
bool QuasilinearDiffusionTerm::CanResonate(
    real_t p,
    real_t xi,
    real_t theta_k,
    int n
) const {
    // Safety check: dispersion and omega_cache must be available
    if (!dispersion || !omega_cache || !dispersion_cached) {
        return true;  // Cannot determine, assume resonance is possible
    }
    
    // Calculate v_parallel and gamma
    real_t gamma = std::sqrt(1.0 + p * p);
    real_t v_par_over_c = (p / gamma) * xi;
    
    // alpha = v_parallel * cos(theta_k)
    real_t alpha = v_par_over_c * Constants::c * std::cos(theta_k);
    
    // beta = 1/gamma
    real_t beta = 1.0 / gamma;
    
    // Get k range from spectrum
    const std::vector<real_t>& k_array = spectrum->getKArray();
    len_t numModes = spectrum->getNumModes();
    
    if (numModes == 0) return false;
    
    real_t k_min = k_array[0];
    real_t k_max = k_array[0];
    for (len_t m = 1; m < numModes; m++) {
        if (k_array[m] < k_min) k_min = k_array[m];
        if (k_array[m] > k_max) k_max = k_array[m];
    }
    
    // Calculate resonant frequency range
    // omega_res = k * alpha + n * omega_ce * beta
    real_t omega_ce = dispersion->getOmegaCE();  // Cyclotron frequency
    real_t omega_res_min = std::min(k_min * alpha + n * omega_ce * beta,
                                    k_max * alpha + n * omega_ce * beta);
    real_t omega_res_max = std::max(k_min * alpha + n * omega_ce * beta,
                                    k_max * alpha + n * omega_ce * beta);
    
    // Get wave spectrum frequency range
    real_t omega_wave_min = omega_cache[0];
    real_t omega_wave_max = omega_cache[0];
    for (len_t m = 1; m < numModes; m++) {
        if (omega_cache[m] < omega_wave_min) omega_wave_min = omega_cache[m];
        if (omega_cache[m] > omega_wave_max) omega_wave_max = omega_cache[m];
    }
    
    // Check for overlap
    // No overlap if: omega_res_max < omega_wave_min OR omega_res_min > omega_wave_max
    if (omega_res_max < omega_wave_min || omega_res_min > omega_wave_max) {
        return false;  // Cannot resonate
    }
    
    return true;  // Possible resonance
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
 * Calculate diffusion coefficients for given radial grid point
 * 
 * Loops over all wave modes and harmonics, finds resonant particles,
 * and accumulates diffusion coefficients.
 */
void QuasilinearDiffusionTerm::calculateDiffusionCoefficients(len_t ir) {
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
    
    // Display progress for this radial point
    std::cerr << "  Assembling diffusion operator at r=" << ir+1 << "/" << grid->GetNr() 
              << " (" << numModes << " modes × " << numHarmonics << " harmonics)..." << std::endl;
    
    // Initialize caches to zero
    for (len_t i = 0; i < (np1 + 1) * np2; i++)
        D_pp_cache[ir * (np1 + 1) * np2 + i] = 0.0;
    
    for (len_t i = 0; i < np1 * (np2 + 1); i++) {
        D_pxi_cache[ir * np1 * (np2 + 1) + i] = 0.0;
        D_xixi_cache[ir * np1 * (np2 + 1) + i] = 0.0;
    }
    
    // Progress tracking
    len_t mode_count = 0;
    
    // Loop over all wave modes
    for (len_t m = 0; m < numModes; m++) {
        real_t k = spectrum->getK(m);
        real_t theta_k = spectrum->getKtheta(m);
        real_t amplitude = spectrum->getAmplitude(m);
        
        // DEBUG: Print amplitude for first few modes
        if (m < 3) {
            std::cerr << "DEBUG calculateDiffusion: mode=" << m << ", amplitude=" << amplitude << std::endl;
        }
        
        if (amplitude < 1e-30) continue;  // Skip negligible modes
        
        mode_count++;
        
        // Display progress every 25% or at completion
        if (mode_count % (numModes / 4) == 0 || m == numModes - 1) {
            int percent = static_cast<int>(mode_count * 100.0 / numModes);
            std::cerr << "    Mode progress: " << percent << "% (" << mode_count << "/" << numModes << ")" << std::endl;
        }
        
        // Use cached omega instead of recalculating
        real_t omega = omega_cache[m];
        if (omega < 0) continue;  // Skip invalid modes
        
        real_t k_parallel = k * std::cos(theta_k);
        
        // Loop over harmonic numbers
        for (int n : harmonicModes) {
            // DEBUG: Print info for first mode and first harmonic
            if (mode_count == 1 && n == 0 && ir == 0) {
                std::cerr << "DEBUG D12: First resonant point - k=" << k << ", theta_k=" << theta_k 
                          << ", k_parallel=" << k_parallel << ", omega=" << omega << std::endl;
            }
            
            // Sample momentum grid points on flux grids
            // D_pp: f1 grid ((np1+1) × np2)
            for (len_t j = 0; j < np2; j++) {
                real_t xi = mg->GetP2(j);  // xi at cell center
                real_t xi_abs = std::abs(xi);
                
                // Avoid singularities at ξ = ±1
                if (xi_abs > 0.999) continue;
                
                for (len_t i = 0; i < np1 + 1; i++) {  // f1 grid: i goes to np1
                    real_t p = mg->GetP1_f(i);  // p at f1 flux face
                    
                    // ⭐ QUADRE-style quick resonance check (fast pre-filter)
                    if (!CanResonate(p, xi, theta_k, n)) {
                        continue;  // Skip this point - cannot resonate with any wave mode
                    }
                    
                    // Use cached omega to find resonant modes (full solver)
                    const std::vector<real_t>& k_array = spectrum->getKArray();
                    std::vector<real_t> resonant_ks = resonanceSolver->findResonantKWithCachedOmega(
                        p, xi, theta_k, n,
                        omega_cache,  // Pre-calculated omega values
                        k_array.data(),  // Pointer to k values array
                        numModes
                    );
                    
                    // Check if current k is near resonance
                    bool is_resonant = false;
                    for (real_t k_res : resonant_ks) {
                        if (std::abs(k - k_res) / k_res < 0.1) {  // Within 10%
                            is_resonant = true;
                            break;
                        }
                    }
                    
                    if (!is_resonant) continue;
                    
                    // DEBUG: Count resonant points
                    static int resonant_count = 0;
                    resonant_count++;
                    if (resonant_count <= 5) {
                        std::cerr << "DEBUG resonant point #" << resonant_count 
                                  << ": mode=" << mode_count << ", n=" << n 
                                  << ", p=" << p << ", xi=" << xi << std::endl;
                    }
                    
                    // Calculate coupling strength |Ψ_{n,k}|²
                    real_t psi_squared = coupling->calculateCouplingStrength(
                        p, xi, k, theta_k, n, *dispersion
                    );
                    
                    // Scale by wave amplitude
                    psi_squared *= amplitude * amplitude;
                    
                    // DEBUG: Check psi_squared after scaling (only for first mode and center grid point)
                    static bool debug_printed = false;
                    if (!debug_printed && mode_count == 0) {
                        std::cerr << "DEBUG psi_squared: mode=" << mode_count << ", amplitude=" << amplitude 
                                  << ", psi_squared=" << psi_squared << std::endl;
                        debug_printed = true;
                    }
                    
                    // Calculate geometric factors
                    real_t gamma = std::sqrt(1.0 + p * p);
                    real_t v_perp_over_c = (p / gamma) * std::sqrt(1.0 - xi * xi);
                    real_t kpar_vperp_over_omega = k_parallel * v_perp_over_c * 3e8 / omega;
                    
                    // Accumulate diffusion coefficients on correct flux grids
                    // D_pp
                    len_t idx_pp = ir * (np1 + 1) * np2 + i * np2 + j;
                    D_pp_cache[idx_pp] += psi_squared;
                }
            }
            
            // D_pξ, D_ξξ: f2 grid (np1 × (np2+1))
            for (len_t j = 0; j < np2 + 1; j++) {  // f2 grid: j goes to np2
                real_t xi = mg->GetP2_f(j);  // xi at f2 flux face
                real_t xi_abs = std::abs(xi);
                
                // Explicitly set ξ = ±1 boundaries to zero
                if (j == 0 || j == np2) {
                    for (len_t i = 0; i < np1; i++) {
                        len_t idx_pxi = ir * np1 * (np2 + 1) + i * (np2 + 1) + j;
                        len_t idx_xixi = ir * np1 * (np2 + 1) + i * (np2 + 1) + j;
                        D_pxi_cache[idx_pxi] = 0.0;
                        D_xixi_cache[idx_xixi] = 0.0;
                    }
                    continue;
                }
                
                // Avoid singularities at ξ = ±1
                if (xi_abs > 0.999) continue;
                
                for (len_t i = 0; i < np1; i++) {  // p at cell center
                    real_t p = mg->GetP1(i);
                    
                    // ⭐ QUADRE-style quick resonance check (fast pre-filter)
                    if (!CanResonate(p, xi, theta_k, n)) {
                        continue;  // Skip this point - cannot resonate with any wave mode
                    }
                    
                    // Use cached omega to find resonant modes (full solver)
                    const std::vector<real_t>& k_array = spectrum->getKArray();
                    std::vector<real_t> resonant_ks = resonanceSolver->findResonantKWithCachedOmega(
                        p, xi, theta_k, n,
                        omega_cache,  // Pre-calculated omega values
                        k_array.data(),  // Pointer to k values array
                        numModes
                    );
                    
                    // Check if current k is near resonance
                    bool is_resonant = false;
                    for (real_t k_res : resonant_ks) {
                        if (std::abs(k - k_res) / k_res < 0.1) {  // Within 10%
                            is_resonant = true;
                            break;
                        }
                    }
                    
                    if (!is_resonant) continue;
                    
                    // DEBUG: Count resonant points
                    static int resonant_count = 0;
                    resonant_count++;
                    if (resonant_count <= 5) {
                        std::cerr << "DEBUG resonant point #" << resonant_count 
                                  << ": mode=" << mode_count << ", n=" << n 
                                  << ", p=" << p << ", xi=" << xi << std::endl;
                    }
                    
                    // Calculate coupling strength |Ψ_{n,k}|²
                    real_t psi_squared = coupling->calculateCouplingStrength(
                        p, xi, k, theta_k, n, *dispersion
                    );
                    
                    // Scale by wave amplitude
                    psi_squared *= amplitude * amplitude;
                    
                    // DEBUG: Check psi_squared after scaling (only for first mode and center grid point)
                    static bool debug_printed = false;
                    if (!debug_printed && mode_count == 0) {
                        std::cerr << "DEBUG psi_squared: mode=" << mode_count << ", amplitude=" << amplitude 
                                  << ", psi_squared=" << psi_squared << std::endl;
                        debug_printed = true;
                    }
                    
                    // Calculate geometric factors
                    real_t gamma = std::sqrt(1.0 + p * p);
                    real_t v = p / gamma;  // Total velocity (normalized to c)
                    real_t v_perp_over_c = v * std::sqrt(1.0 - xi * xi);
                    
                    // CORRECTED geometric factors according to Guo 2018 theory:
                    // D_pξ ∝ -√(1-ξ²) · (ξ - k_∥v/ω)
                    // D_ξξ ∝ (ξ - k_∥v/ω)²
                    real_t xi_minus_kpar_v_over_omega = xi - k_parallel * v * 3e8 / omega;
                    
                    // DEBUG: Check geometric factor
                    if (mode_count == 1 && n == 0 && ir == 0 && i == np1/2 && j == np2/2) {
                        std::cerr << "DEBUG D12 geometry: xi=" << xi << ", k_par*v/omega=" << k_parallel * v * 3e8 / omega
                                  << ", (xi - k_par*v/omega)=" << xi_minus_kpar_v_over_omega
                                  << ", sqrt(1-xi^2)=" << std::sqrt(1.0 - xi*xi) << std::endl;
                    }
                    
                    // Accumulate diffusion coefficients
                    // D_pξ (with negative sign from theory)
                    len_t idx_pxi = ir * np1 * (np2 + 1) + i * (np2 + 1) + j;
                    real_t contribution_pxi = -psi_squared * std::sqrt(1.0 - xi * xi) * xi_minus_kpar_v_over_omega;
                    D_pxi_cache[idx_pxi] += contribution_pxi;
                    
                    // DEBUG: Track D12 contributions
                    if (mode_count <= 5 && std::abs(contribution_pxi) > 1e-30) {
                        std::cerr << "DEBUG D12 contribution: mode=" << mode_count << ", n=" << n 
                                  << ", p=" << p << ", xi=" << xi 
                                  << ", psi^2=" << psi_squared 
                                  << ", geom_factor=" << xi_minus_kpar_v_over_omega
                                  << ", sqrt(1-xi^2)=" << std::sqrt(1.0 - xi*xi)
                                  << ", contrib=" << contribution_pxi << std::endl;
                    }
                    
                    // Count total D12 contributions for debugging
                    static int d12_contribution_count = 0;
                    if (std::abs(contribution_pxi) > 1e-30) {
                        d12_contribution_count++;
                        if (d12_contribution_count <= 10) {
                            std::cerr << "D12_NONZERO: mode=" << mode_count << ", n=" << n 
                                      << ", ir=" << ir << ", i=" << i << ", j=" << j
                                      << ", val=" << contribution_pxi << std::endl;
                        }
                    }
                    
                    // D_ξξ
                    len_t idx_xixi = ir * np1 * (np2 + 1) + i * (np2 + 1) + j;
                    D_xixi_cache[idx_xixi] += psi_squared * xi_minus_kpar_v_over_omega * xi_minus_kpar_v_over_omega;
                }
            }
        }
    }
    
    // Save current amplitudes for change detection
    if (!last_amplitudes || num_modes_tracked != numModes) {
        if (last_amplitudes) delete[] last_amplitudes;
        last_amplitudes = new real_t[numModes];
        num_modes_tracked = numModes;
    }
    for (len_t m = 0; m < numModes; m++) {
        last_amplitudes[m] = spectrum->getAmplitude(m);
    }
    
    // DEBUG: Check cache values
    static bool debug_cache_checked = false;
    if (!debug_cache_checked && ir == 0) {
        real_t max_D_pp = 0.0, sum_D_pp = 0.0;
        len_t nonzero_D_pp = 0;
        for (len_t i = 0; i < (np1 + 1) * np2; i++) {
            real_t val = std::abs(D_pp_cache[ir * (np1 + 1) * np2 + i]);
            if (val > max_D_pp) max_D_pp = val;
            sum_D_pp += val;
            if (val > 1e-30) nonzero_D_pp++;
        }
        std::cerr << "DEBUG cache check: D_pp max=" << max_D_pp << ", sum=" << sum_D_pp 
                  << ", nonzero=" << nonzero_D_pp << "/" << ((np1+1)*np2) << std::endl;
        debug_cache_checked = true;
    }
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
void QuasilinearDiffusionTerm::Rebuild(const real_t /*t*/, const real_t /*dt*/, FVM::UnknownQuantityHandler* /*uqh*/) {
    // Check if wave spectrum amplitude has changed
    bool amplitude_changed = false;
    if (spectrum && operator_cached) {
        len_t numModes = spectrum->getNumK() * spectrum->getNumKtheta();
        for (len_t m = 0; m < numModes; m++) {
            real_t current_amp = spectrum->getAmplitude(m);
            if (std::abs(current_amp - last_amplitudes[m]) > 1e-30) {
                amplitude_changed = true;
                break;
            }
        }
    }
    
    // If amplitude changed, invalidate cache and recalculate
    if (amplitude_changed) {
        std::cerr << "Wave amplitude changed, recalculating diffusion coefficients..." << std::endl;
        operator_cached = false;
    }
    
    const len_t nr = grid->GetNr();
    
    // Allocate or reallocate cache if needed (first time or grid size changed)
    if (!operator_cached || nr != nr_cached) {
        if (D_pp_cache) delete[] D_pp_cache;
        if (D_pxi_cache) delete[] D_pxi_cache;
        if (D_xixi_cache) delete[] D_xixi_cache;
        
        auto *mg = grid->GetMomentumGrid(0);
        np1_cached = mg->GetNp1();
        np2_cached = mg->GetNp2();
        
        D_pp_cache = new real_t[nr * (np1_cached + 1) * np2_cached];
        D_pxi_cache = new real_t[nr * np1_cached * (np2_cached + 1)];
        D_xixi_cache = new real_t[nr * np1_cached * (np2_cached + 1)];
        
        nr_cached = nr;
    }
    
    // Calculate or load diffusion coefficients for each radial point (only if not cached)
    if (!operator_cached) {
        if (use_precomputed_matrix && matrix_loader) {
        // Use pre-computed matrix with amplitude scaling
        for (len_t ir = 0; ir < nr; ir++) {
            auto *mg = grid->GetMomentumGrid(ir);
            const len_t np1 = mg->GetNp1();
            const len_t np2 = mg->GetNp2();
            
            // Zero out caches
            memset(D_pp_cache + ir * (np1 + 1) * np2, 0, sizeof(real_t) * (np1 + 1) * np2);
            memset(D_pxi_cache + ir * np1 * (np2 + 1), 0, sizeof(real_t) * np1 * (np2 + 1));
            memset(D_xixi_cache + ir * np1 * (np2 + 1), 0, sizeof(real_t) * np1 * (np2 + 1));
            
            // Fill from pre-computed matrix (already scaled by amplitude)
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1 + 1; i++) {
                    if (i < np1) {  // Only access valid indices
                        D_pp_cache[ir * (np1 + 1) * np2 + i * np2 + j] = 
                            matrix_loader->getDppScaled(i, j);
                    }
                }
            }
            
            for (len_t j = 0; j < np2 + 1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    D_pxi_cache[ir * np1 * (np2 + 1) + i * (np2 + 1) + j] = 
                        matrix_loader->getDpxiScaled(i, j);
                    D_xixi_cache[ir * np1 * (np2 + 1) + i * (np2 + 1) + j] = 
                        matrix_loader->getDxixiScaled(i, j);
                }
            }
        }
        
        std::cerr << "✓ Quasi-linear diffusion operator loaded from pre-computed matrix" << std::endl;
        
        // ====================================================================
        // NaN Check: Validate all loaded diffusion coefficients
        // ====================================================================
        bool has_nan = false;
        len_t nan_count_pp = 0, nan_count_pxi = 0, nan_count_xixi = 0;
        
        for (len_t ir = 0; ir < nr; ir++) {
            auto *mg = grid->GetMomentumGrid(ir);
            const len_t np1 = mg->GetNp1();
            const len_t np2 = mg->GetNp2();
            
            // Check D_pp
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1 + 1; i++) {
                    real_t val = D_pp_cache[ir * (np1 + 1) * np2 + i * np2 + j];
                    if (std::isnan(val) || std::isinf(val)) {
                        has_nan = true;
                        nan_count_pp++;
                        if (nan_count_pp <= 5) {
                            std::cerr << "ERROR: NaN/Inf in D_pp at ir=" << ir 
                                     << ", i=" << i << ", j=" << j 
                                     << ", value=" << val << std::endl;
                        }
                    }
                }
            }
            
            // Check D_pξ and D_ξξ
            for (len_t j = 0; j < np2 + 1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    real_t val_pxi = D_pxi_cache[ir * np1 * (np2 + 1) + i * (np2 + 1) + j];
                    real_t val_xixi = D_xixi_cache[ir * np1 * (np2 + 1) + i * (np2 + 1) + j];
                    
                    if (std::isnan(val_pxi) || std::isinf(val_pxi)) {
                        has_nan = true;
                        nan_count_pxi++;
                        if (nan_count_pxi <= 5) {
                            std::cerr << "ERROR: NaN/Inf in D_pξ at ir=" << ir 
                                     << ", i=" << i << ", j=" << j 
                                     << ", value=" << val_pxi << std::endl;
                        }
                    }
                    
                    if (std::isnan(val_xixi) || std::isinf(val_xixi)) {
                        has_nan = true;
                        nan_count_xixi++;
                        if (nan_count_xixi <= 5) {
                            std::cerr << "ERROR: NaN/Inf in D_ξξ at ir=" << ir 
                                     << ", i=" << i << ", j=" << j 
                                     << ", value=" << val_xixi << std::endl;
                        }
                    }
                }
            }
        }
        
        if (has_nan) {
            std::cerr << "\n========================================" << std::endl;
            std::cerr << "FATAL ERROR: Quasi-linear diffusion coefficients contain NaN/Inf!" << std::endl;
            std::cerr << "  D_pp NaN/Inf count:   " << nan_count_pp << std::endl;
            std::cerr << "  D_pξ NaN/Inf count:   " << nan_count_pxi << std::endl;
            std::cerr << "  D_ξξ NaN/Inf count:   " << nan_count_xixi << std::endl;
            std::cerr << "\nPossible causes:" << std::endl;
            std::cerr << "  1. Pre-computed HDF5 file contains invalid data" << std::endl;
            std::cerr << "  2. Wave amplitude scaling produced overflow" << std::endl;
            std::cerr << "  3. Grid mismatch between HDF5 and DREAM" << std::endl;
            std::cerr << "  4. Numerical issues in Python precomputation" << std::endl;
            std::cerr << "\nAborting simulation to prevent numerical instability." << std::endl;
            std::cerr << "========================================\n" << std::endl;
            throw DREAM::DREAMException("Quasi-linear diffusion coefficients contain NaN/Inf values");
        }
        
        std::cerr << "✓ All diffusion coefficients validated (no NaN/Inf detected)" << std::endl;
        } else {
            // Use REVERSE resonance solving method for grid-independent results
            for (len_t ir = 0; ir < nr; ir++) {
                calculateDiffusionCoefficientsReverse(ir);
            }
            
            std::cerr << "✓ Quasi-linear diffusion operator computed on-the-fly (REVERSE method)" << std::endl;
        }
    }
    
    // Mark operator as cached - no need to recalculate in future time steps
    operator_cached = true;
    
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
                real_t D_pp_val = getD_pp(ir, i, j);
                
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
                real_t D_xixi_val = getD_xixi(ir, i, j);
                
                // Final safety check
                if (std::isnan(D_xixi_val) || std::isinf(D_xixi_val)) {
                    std::cerr << "FATAL: Attempting to fill NaN/Inf into D22 at ir=" << ir 
                             << ", i=" << i << ", j=" << j << std::endl;
                    throw DREAM::DREAMException("NaN/Inf detected in D_xixi during FVM matrix assembly");
                }
                
                D22(ir, i, j) += D_xixi_val;
            }
        }
        
        // D12 = D21 = D_pξ (cross diffusion from quasi-linear theory)
        // Fill from pre-computed cache (Zehua Guo Eq. 11b)
        for (len_t j = 0; j < np2 + 1; j++) {
            for (len_t i = 0; i < np1; i++) {
                real_t D_pxi_val = getD_pxi(ir, i, j);
                
                // Final safety check for cross terms
                if (std::isnan(D_pxi_val) || std::isinf(D_pxi_val)) {
                    std::cerr << "FATAL: Attempting to fill NaN/Inf into D12/D21 at ir=" << ir 
                             << ", i=" << i << ", j=" << j << std::endl;
                    throw DREAM::DREAMException("NaN/Inf detected in D_pxi during FVM matrix assembly");
                }
                
                D12(ir, i, j) += D_pxi_val;
                D21(ir, i, j) += D_pxi_val;  // Symmetric tensor
            }
        }
        
        // ===== VALIDATION: Check if diffusion coefficients are non-zero (only on first rebuild, first radial point) =====
        // NOTE: This must be done AFTER filling the FVM matrices
        if (!first_rebuild_done && ir == 0) {
            auto *mg = grid->GetMomentumGrid(ir);
            const len_t np1 = mg->GetNp1();
            const len_t np2 = mg->GetNp2();
            
            real_t max_D11 = 0.0, max_D22 = 0.0, max_D12 = 0.0;
            real_t sum_D11 = 0.0, sum_D22 = 0.0, sum_D12 = 0.0;
            len_t count_nonzero_D11 = 0, count_nonzero_D22 = 0, count_nonzero_D12 = 0;
            
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1 + 1; i++) {
                    real_t val = std::abs(D11(ir, i, j));
                    if (val > max_D11) max_D11 = val;
                    sum_D11 += val;
                    if (val > 1e-30) count_nonzero_D11++;
                }
            }
            
            for (len_t j = 0; j < np2 + 1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    real_t val = std::abs(D22(ir, i, j));
                    if (val > max_D22) max_D22 = val;
                    sum_D22 += val;
                    if (val > 1e-30) count_nonzero_D22++;
                    
                    if (j < np2) {
                        val = std::abs(D12(ir, i, j));
                        if (val > max_D12) max_D12 = val;
                        sum_D12 += val;
                        if (val > 1e-30) count_nonzero_D12++;
                    }
                }
            }
            
            std::cerr << "\n=== QUASILINEAR DIFFUSION VALIDATION (r=1/" << grid->GetNr() << ") ===" << std::endl;
            std::cerr << "D11 (p-p): max=" << max_D11 << ", sum=" << sum_D11 
                      << ", nonzero=" << count_nonzero_D11 << "/" << ((np1+1)*np2) << std::endl;
            std::cerr << "D22 (ξ-ξ): max=" << max_D22 << ", sum=" << sum_D22 
                      << ", nonzero=" << count_nonzero_D22 << "/" << (np1*(np2+1)) << std::endl;
            std::cerr << "D12 (p-ξ): max=" << max_D12 << ", sum=" << sum_D12 
                      << ", nonzero=" << count_nonzero_D12 << "/" << (np1*np2) << std::endl;
            
            if (max_D11 < 1e-30 && max_D22 < 1e-30 && max_D12 < 1e-30) {
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
    
    // Display progress
    std::cerr << "  Assembling diffusion operator (REVERSE method) at r=" << ir+1 << "/" << grid->GetNr() 
              << " (" << numModes << " modes × " << numHarmonics << " harmonics)..." << std::endl;
    
    // Initialize caches to zero
    for (len_t i = 0; i < (np1 + 1) * np2; i++)
        D_pp_cache[ir * (np1 + 1) * np2 + i] = 0.0;
    
    for (len_t i = 0; i < np1 * (np2 + 1); i++) {
        D_pxi_cache[ir * np1 * (np2 + 1) + i] = 0.0;
        D_xixi_cache[ir * np1 * (np2 + 1) + i] = 0.0;
    }
    
    // Resonance width parameter (fraction of grid spacing)
    // This controls how wide the resonance peak is spread
    const real_t resonance_width_factor = 0.5;  // 50% of grid spacing
    
    // Progress tracking
    len_t mode_count = 0;
    int total_resonance_points = 0;
    
    // Loop over all wave modes
    for (len_t m = 0; m < numModes; m++) {
        real_t k = spectrum->getK(m);
        real_t theta_k = spectrum->getKtheta(m);
        real_t amplitude = spectrum->getAmplitude(m);
        
        if (amplitude < 1e-30) continue;  // Skip negligible modes
        
        mode_count++;
        
        // Display progress every 25%
        if (mode_count % (numModes / 4) == 0 || m == numModes - 1) {
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
                
                // DEBUG: Check if any resonance found
                static std::map<int, int> debug_count_per_harmonic;
                if (debug_count_per_harmonic.find(n) == debug_count_per_harmonic.end()) {
                    debug_count_per_harmonic[n] = 0;
                }
                
                // Only print for xi close to ±1 (within 0.1)
                bool xi_near_boundary = (std::abs(std::abs(xi) - 1.0) < 0.1);
                
                if (xi_near_boundary && debug_count_per_harmonic[n] < 3 && ir == 0 && m == 0) {
                    std::cerr << "DEBUG findResonantP: mode=" << m << ", k=" << k << ", theta_k=" << theta_k 
                              << ", n=" << n << ", xi=" << xi 
                              << ", found " << resonant_ps.size() << " solutions" << std::endl;
                    for (auto p_res : resonant_ps) {
                        std::cerr << "  -> p_res=" << p_res << std::endl;
                    }
                    debug_count_per_harmonic[n]++;
                }
                
                if (resonant_ps.empty()) continue;  // No resonance for this (xi, n)
                
                // Process each resonant momentum
                for (real_t p_res : resonant_ps) {
                    total_resonance_points++;
                    
                    // DEBUG: Print first few resonant points
                    if (total_resonance_points <= 3 && ir == 0) {
                        // std::cerr << "DEBUG reverse resonance: mode=" << m << ", n=" << n 
                        //           << ", xi=" << xi << ", p_res=" << p_res << std::endl;
                    }
                    
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
                    
                    // Calculate coupling strength at resonant point
                    real_t psi_squared = coupling->calculateCouplingStrength(
                        p_res, xi, k, theta_k, n, *dispersion
                    );
                    
                    // Scale by wave amplitude
                    psi_squared *= amplitude * amplitude;
                    
                    // Apply resonance weight
                    psi_squared *= weight_sum;  // Total weight for normalization
                    
                    // Calculate geometric factors at resonant point
                    real_t gamma = std::sqrt(1.0 + p_res * p_res);
                    real_t v_par_over_c = (p_res / gamma) * xi;
                    real_t c_val = 2.99792458e8;  // Speed of light (m/s)
                    real_t xi_minus_kpar_v_over_omega = xi - k_parallel * v_par_over_c * c_val / omega;
                    
                    // Accumulate to D_pp (at f1 grid points)
                    // Distribute to i_low and i_high with weights
                    len_t idx_pp_low = ir * (np1 + 1) * np2 + i_low * np2 + j;
                    len_t idx_pp_high = ir * (np1 + 1) * np2 + i_high * np2 + j;
                    
                    D_pp_cache[idx_pp_low] += psi_squared * weight_low;
                    D_pp_cache[idx_pp_high] += psi_squared * weight_high;
                    
                    // Accumulate to D_pξ and D_ξξ (at f2 grid points)
                    // For f2 grid, p index ranges from 0 to np1-1
                    // Need to clamp i_low and i_high to valid range
                    len_t i_low_f2 = std::min(i_low, np1 - 1);
                    len_t i_high_f2 = std::min(i_high, np1 - 1);
                    
                    // BOUNDS CHECK: Ensure indices are valid
                    if (i_low_f2 >= np1 || j >= (np2 + 1)) {
                        std::cerr << "ERROR: Invalid index for D_pxi/D_xixi: i_low_f2=" << i_low_f2 
                                  << ", j=" << j << ", np1=" << np1 << ", np2=" << np2 << std::endl;
                        throw DREAMException("Array index out of bounds in calculateDiffusionCoefficientsReverse");
                    }
                    if (i_high_f2 >= np1 || j >= (np2 + 1)) {
                        std::cerr << "ERROR: Invalid index for D_pxi/D_xixi: i_high_f2=" << i_high_f2 
                                  << ", j=" << j << ", np1=" << np1 << ", np2=" << np2 << std::endl;
                        throw DREAMException("Array index out of bounds in calculateDiffusionCoefficientsReverse");
                    }
                    
                    len_t idx_pxi_low = ir * np1 * (np2 + 1) + i_low_f2 * (np2 + 1) + j;
                    len_t idx_pxi_high = ir * np1 * (np2 + 1) + i_high_f2 * (np2 + 1) + j;
                    len_t idx_xixi_low = ir * np1 * (np2 + 1) + i_low_f2 * (np2 + 1) + j;
                    len_t idx_xixi_high = ir * np1 * (np2 + 1) + i_high_f2 * (np2 + 1) + j;
                    
                    real_t contrib_pxi = -psi_squared * std::sqrt(1.0 - xi * xi) * xi_minus_kpar_v_over_omega;
                    real_t contrib_xixi = psi_squared * xi_minus_kpar_v_over_omega * xi_minus_kpar_v_over_omega;
                    
                    D_pxi_cache[idx_pxi_low] += contrib_pxi * weight_low;
                    D_pxi_cache[idx_pxi_high] += contrib_pxi * weight_high;
                    
                    D_xixi_cache[idx_xixi_low] += contrib_xixi * weight_low;
                    D_xixi_cache[idx_xixi_high] += contrib_xixi * weight_high;
                }
            }
        }
    }
    
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
