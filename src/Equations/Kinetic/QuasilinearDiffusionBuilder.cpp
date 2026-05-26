/**
 * Builder for quasilinear diffusion equation terms.
 * 
 * This file contains the logic for constructing QuasilinearDiffusionTerm objects
 * from Settings configuration, supporting both pre-computed matrix and on-the-fly
 * computation modes.
 */

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "DREAM/Equations/Kinetic/QuasilinearDiffusionBuilder.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/WaveSpectrum.hpp"
#include "DREAM/Equations/Kinetic/WhistlerDispersion.hpp"
#include "DREAM/Equations/Kinetic/ResonanceSolver.hpp"
#include "DREAM/Equations/Kinetic/WaveParticleCoupling.hpp"
#include "DREAM/Equations/Kinetic/QuasilinearDiffusionTerm.hpp"
#include "DREAM/DREAMException.hpp"

using namespace DREAM;

/**
 * Construct a quasilinear diffusion term from settings.
 * 
 * Supports two modes:
 * 1. Pre-computed matrix mode: Loads diffusion coefficients from HDF5 file
 * 2. On-the-fly computation mode: Calculates coefficients using wave-particle interaction theory
 */
QuasilinearDiffusionTerm *DREAM::ConstructQuasilinearDiffusionTerm(
    Settings *s, 
    const std::string& mod, 
    FVM::Grid *grid
) {
    std::cerr << "DEBUG: ConstructQuasilinearDiffusionTerm called with mod=" << mod << std::endl;
    
    enum OptionConstants::ql_diffusion_mode ql_mode =
        (enum OptionConstants::ql_diffusion_mode)s->GetInteger(mod + "/quasilinearmode");
    
    std::cerr << "DEBUG: ql_mode = " << ql_mode << std::endl;
    
    if (ql_mode == OptionConstants::QL_DIFFUSION_MODE_NEGLECT)
        return nullptr;
    
    // Check if using pre-computed matrix
    bool use_precomputed = s->GetInteger(mod + "/quasilinear/use_precomputed_matrix") != 0;
    
    if (use_precomputed) {
        // ====== Mode A: Pre-computed matrix from HDF5 ======
        
        // Safely get precomputed file path
        std::string hdf5_file;
        try {
            hdf5_file = s->GetString(mod + "/quasilinear/precomputed_file");
        } catch (...) {
            throw SettingsException(
                "%s/quasilinear/precomputed_file is required when use_precomputed_matrix=1",
                mod.c_str()
            );
        }
        
        if (hdf5_file.empty()) {
            throw SettingsException(
                "%s/quasilinear/precomputed_file must not be empty when use_precomputed_matrix=1",
                mod.c_str()
            );
        }
        
        real_t amplitude = s->GetReal(mod + "/quasilinear/amplitude");
        real_t start_inject_time = s->GetReal(mod + "/quasilinear/start_inject_time");
        real_t inject_cycle_duration = s->GetReal(mod + "/quasilinear/inject_cycle_duration");
        
        QuasilinearDiffusionTerm *qlTerm = new QuasilinearDiffusionTerm(
            grid, hdf5_file, amplitude, start_inject_time, inject_cycle_duration
        );
        
        std::cerr << "✓ Quasilinear diffusion enabled with pre-computed matrix: " 
                  << hdf5_file << std::endl;
        std::cerr << "  Initial amplitude: " << amplitude << std::endl;
        
        return qlTerm;
        
    } else {
        // ====== Mode B: On-the-fly computation ======
        
        // Get plasma parameters for dispersion relation
        real_t B0 = s->GetReal("radialgrid/B0");  // Magnetic field strength from radial grid
        std::cerr << "DEBUG: Read B0 from settings: B0=" << B0 << " T" << std::endl;
        
        // Get electron density - try to read from first ion species
        real_t n_e = 0.0;
        bool density_found = false;
        
        // Try common ion density paths
        const char* ion_density_paths[] = {
            "eqsys/n_i/ion/0/n",
            "eqsys/n_i/0/n"
        };
        
        for (const auto& path : ion_density_paths) {
            try {
                real_t n_i = s->GetReal(path);
                real_t Z = 1.0;  // Default charge state
                
                // Try to get Z from the same ion
                std::string z_path = std::string(path).substr(0, std::string(path).rfind('/')) + "/Z";
                try {
                    Z = s->GetReal(z_path.c_str());
                } catch (...) {
                    // Use default Z=1
                }
                
                n_e = Z * n_i;
                density_found = true;
                std::cerr << "DEBUG: Read ion density from " << path << ": n_i=" << n_i << ", Z=" << Z << ", n_e=" << n_e << std::endl;
                break;
            } catch (...) {
                // Try next path
            }
        }
        
        if (!density_found || n_e <= 0) {
            std::cerr << "Warning: Could not read ion density, using default n_e = 5e18 m^-3" << std::endl;
            n_e = 5e18;  // More reasonable default based on typical simulations
        }
        
        // Create wave spectrum
        len_t num_k = s->GetInteger(mod + "/quasilinear/num_k");
        len_t num_ktheta = s->GetInteger(mod + "/quasilinear/num_ktheta");
        
        WaveSpectrum *spectrum = new WaveSpectrum(num_k, num_ktheta);
        
        // Set spectrum type and parameters
        enum OptionConstants::wave_spectrum_type spec_type =
            (enum OptionConstants::wave_spectrum_type)s->GetInteger(mod + "/quasilinear/spectrum_type");
        
        if (spec_type == OptionConstants::WAVE_SPECTRUM_UNIFORM) {
            real_t k_min = s->GetReal(mod + "/quasilinear/k_min");
            real_t k_max = s->GetReal(mod + "/quasilinear/k_max");
            real_t ktheta_min = s->GetReal(mod + "/quasilinear/ktheta_min");
            real_t ktheta_max = s->GetReal(mod + "/quasilinear/ktheta_max");
            
            spectrum->setUniformSpectrum(k_min, k_max, ktheta_min, ktheta_max);
        } else {
            std::cerr << "Warning: Only uniform spectrum type is currently supported for quasilinear diffusion." << std::endl;
            spectrum->setUniformSpectrum(35.0, 45.0, 0.1, 0.3);  // Default values
        }
        
        // Set amplitude for all modes
        real_t amplitude = s->GetReal(mod + "/quasilinear/amplitude");
        real_t start_inject_time = s->GetReal(mod + "/quasilinear/start_inject_time");
        real_t inject_cycle_duration = s->GetReal(mod + "/quasilinear/inject_cycle_duration");
        std::cerr << "DEBUG: Read amplitude from settings: amplitude=" << amplitude << std::endl;
        len_t num_modes = spectrum->getNumModes();
        for (len_t m = 0; m < num_modes; m++) {
            spectrum->setAmplitude(m, amplitude);
        }
        std::cerr << "DEBUG: Set amplitude for " << num_modes << " modes" << std::endl;
        
        // Create dispersion relation solver
        bool use_simple_dispersion = s->GetInteger(mod + "/quasilinear/use_simple_dispersion") != 0;
        std::cerr << "DEBUG: Creating WhistlerDispersion with B0=" << B0 << ", n_e=" << n_e << std::endl;
        WhistlerDispersion *dispersion = new WhistlerDispersion(B0, n_e, 1836.0, 1.0, use_simple_dispersion);
        
        if (use_simple_dispersion) {
            std::cerr << "  Using simplified whistler dispersion relation (ω = k|k_∥| * w)" << std::endl;
        } else {
            std::cerr << "  Using full PDRF dispersion relation solver" << std::endl;
        }
        std::cerr << "DEBUG: WhistlerDispersion created successfully" << std::endl;
        
        // Calculate plasma frequencies
        constexpr real_t e_charge = 1.60217662e-19;   // C
        constexpr real_t m_electron = 9.1094e-31;     // kg
        constexpr real_t epsilon_0 = 8.854187817e-12; // F/m
        
        real_t omega_ce = e_charge * B0 / m_electron;
        real_t omega_pe = std::sqrt(n_e * e_charge * e_charge / (m_electron * epsilon_0));
        
        // Create resonance solver
        ResonanceSolver *resonanceSolver = new ResonanceSolver(omega_ce);
        
        // Create wave-particle coupling calculator
        WaveParticleCoupling *coupling = new WaveParticleCoupling(omega_pe, omega_ce, n_e);
        
        // Determine harmonic modes to include
        enum OptionConstants::ql_harmonic_mode hmode =
            (enum OptionConstants::ql_harmonic_mode)s->GetInteger(mod + "/quasilinear/harmonic_mode");
        
        std::vector<int> harmonicModes;
        if (hmode == OptionConstants::QL_HARMONIC_N_MINUS_1) {
            harmonicModes.push_back(-1);
        } else if (hmode == OptionConstants::QL_HARMONIC_N_PLUS_1) {
            harmonicModes.push_back(+1);
        } else if (hmode == OptionConstants::QL_HARMONIC_BOTH) {
            harmonicModes.push_back(-2);
            harmonicModes.push_back(-1);
            harmonicModes.push_back(0);
            harmonicModes.push_back(+1);
            harmonicModes.push_back(+2);
        }
        
        // Create quasilinear diffusion term
        QuasilinearDiffusionTerm *qlTerm = new QuasilinearDiffusionTerm(
            grid, spectrum, dispersion, resonanceSolver, coupling, harmonicModes,
            start_inject_time, inject_cycle_duration
        );
        
        std::cerr << "Quasilinear diffusion enabled with " << num_modes << " wave modes." << std::endl;
        std::cerr << "  Harmonic modes: n = ";
        for (size_t i = 0; i < harmonicModes.size(); i++) {
            std::cerr << harmonicModes[i];
            if (i < harmonicModes.size() - 1) std::cerr << ", ";
        }
        std::cerr << std::endl;
        
        return qlTerm;
    }
}
