/**
 * Builder for quasilinear diffusion equation terms.
 *
 * Constructs QuasilinearDiffusionTerm from Settings, using on-the-fly
 * resonance-based computation of diffusion coefficients.
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

QuasilinearDiffusionTerm *DREAM::ConstructQuasilinearDiffusionTerm(
    Settings *s, 
    const std::string& mod, 
    FVM::Grid *grid
) {
    enum OptionConstants::ql_diffusion_mode ql_mode =
        (enum OptionConstants::ql_diffusion_mode)s->GetInteger(mod + "/quasilinearmode");
    
    if (ql_mode == OptionConstants::QL_DIFFUSION_MODE_NEGLECT)
        return nullptr;
    
    // Read plasma parameters
    real_t B0 = s->GetReal("radialgrid/B0");
    real_t n_e = s->GetReal(mod + "/quasilinear/density");
    
    // Create wave spectrum
    len_t num_k = s->GetInteger(mod + "/quasilinear/num_k");
    len_t num_ktheta = s->GetInteger(mod + "/quasilinear/num_ktheta");
    
    WaveSpectrum *spectrum = new WaveSpectrum(num_k, num_ktheta);
    
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
        spectrum->setUniformSpectrum(35.0, 45.0, 0.1, 0.3);
    }
    
    // Set amplitude for all modes
    real_t amplitude = s->GetReal(mod + "/quasilinear/amplitude");
    real_t start_inject_time = s->GetReal(mod + "/quasilinear/start_inject_time");
    real_t inject_cycle_duration = s->GetReal(mod + "/quasilinear/inject_cycle_duration");
    real_t ramp_time = s->GetReal(mod + "/quasilinear/ramp_time");
    len_t num_modes = spectrum->getNumModes();
    for (len_t m = 0; m < num_modes; m++) {
        spectrum->setAmplitude(m, amplitude);
    }
    
    // Create dispersion relation solver
    bool use_simple_dispersion = s->GetInteger(mod + "/quasilinear/use_simple_dispersion") != 0;
    WhistlerDispersion *dispersion = new WhistlerDispersion(B0, n_e, 1836.0, 1.0, use_simple_dispersion);
    
    if (use_simple_dispersion) {
        std::cerr << "  Using simplified whistler dispersion relation" << std::endl;
    } else {
        std::cerr << "  Using full stix dispersion relation solver" << std::endl;
    }
    
    // Calculate plasma frequencies
    constexpr real_t e_charge = 1.60217662e-19;
    constexpr real_t m_electron = 9.1094e-31;
    constexpr real_t epsilon_0 = 8.854187817e-12;
    
    real_t omega_ce = e_charge * B0 / m_electron;
    real_t omega_pe = std::sqrt(n_e * e_charge * e_charge / (m_electron * epsilon_0));
    
    // Create resonance solver and wave-particle coupling
    ResonanceSolver *resonanceSolver = new ResonanceSolver(omega_ce);
    WaveParticleCoupling *coupling = new WaveParticleCoupling(omega_pe, omega_ce, n_e);
    
    // Determine harmonic modes
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
    
    QuasilinearDiffusionTerm *qlTerm = new QuasilinearDiffusionTerm(
        grid, spectrum, dispersion, resonanceSolver, coupling, harmonicModes,
        start_inject_time, inject_cycle_duration, ramp_time
    );
    
    std::cerr << "Quasilinear diffusion enabled: " << num_modes << " modes, "
              << "n_e=" << n_e << " m^-3, B0=" << B0 << " T" << std::endl;
    std::cerr << "  Harmonic modes: n = ";
    for (size_t i = 0; i < harmonicModes.size(); i++) {
        std::cerr << harmonicModes[i];
        if (i < harmonicModes.size() - 1) std::cerr << ", ";
    }
    std::cerr << std::endl;
    
    return qlTerm;
}
