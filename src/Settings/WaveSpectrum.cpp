/**
 * Implementation of the WaveSpectrum class for managing
 * wave spectrum parameters in quasilinear diffusion.
 */

#include "DREAM/Settings/WaveSpectrum.hpp"
#include "DREAM/Constants.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace DREAM;

/**
 * Constructor
 */
WaveSpectrum::WaveSpectrum(len_t nk, len_t nkt)
    : num_k(nk), num_ktheta(nkt), num_modes(nk * nkt),
      spectrum_type(OptionConstants::WAVE_SPECTRUM_UNIFORM) {
    
    initializeArrays();
}

/**
 * Destructor
 */
WaveSpectrum::~WaveSpectrum() {
    // Vectors will be automatically cleaned up
}

/**
 * Initialize arrays with given sizes
 */
void WaveSpectrum::initializeArrays() {
    k_array.resize(num_modes, 0.0);
    ktheta_array.resize(num_modes, 0.0);
    amplitude.resize(num_modes, 0.0);
}

/**
 * Set uniform spectrum in (k, theta_k) space
 */
void WaveSpectrum::setUniformSpectrum(real_t k_min, real_t k_max, 
                                     real_t ktheta_min, real_t ktheta_max) {
    spectrum_type = OptionConstants::WAVE_SPECTRUM_UNIFORM;
    
    uniform_params.k_min = k_min;
    uniform_params.k_max = k_max;
    uniform_params.ktheta_min = ktheta_min;
    uniform_params.ktheta_max = ktheta_max;
    
    // Generate uniform grid
    real_t dk = (k_max - k_min) / (num_k - 1);
    real_t dktheta = (ktheta_max - ktheta_min) / (num_ktheta - 1);
    
    for (len_t i = 0; i < num_ktheta; i++) {
        for (len_t j = 0; j < num_k; j++) {
            len_t idx = i * num_k + j;
            k_array[idx] = k_min + j * dk;
            ktheta_array[idx] = ktheta_min + i * dktheta;
        }
    }
    
    std::cout << "WaveSpectrum: Set uniform spectrum with " 
              << num_k << " x " << num_ktheta << " = " << num_modes << " modes" << std::endl;
    std::cout << "  k range: [" << k_min << ", " << k_max << "]" << std::endl;
    std::cout << "  theta_k range: [" << ktheta_min << ", " << ktheta_max << "]" << std::endl;
}

/**
 * Set Gaussian spectrum centered at (k0, ktheta0)
 */
void WaveSpectrum::setGaussianSpectrum(real_t k0, real_t ktheta0, 
                                      real_t delta_k, real_t delta_ktheta) {
    spectrum_type = OptionConstants::WAVE_SPECTRUM_GAUSSIAN;
    
    gaussian_params.k0 = k0;
    gaussian_params.ktheta0 = ktheta0;
    gaussian_params.delta_k = delta_k;
    gaussian_params.delta_ktheta = delta_ktheta;
    
    // Generate grid covering ±3 sigma
    real_t k_min = k0 - 3 * delta_k;
    real_t k_max = k0 + 3 * delta_k;
    real_t ktheta_min = ktheta0 - 3 * delta_ktheta;
    real_t ktheta_max = ktheta0 + 3 * delta_ktheta;
    
    // Ensure positive values
    k_min = std::max(k_min, 1e-10);
    ktheta_min = std::max(ktheta_min, 1e-10);
    
    setUniformSpectrum(k_min, k_max, ktheta_min, ktheta_max);
    
    // Set Gaussian amplitudes
    real_t norm = 1.0 / (std::sqrt(2*M_PI) * delta_k * std::sqrt(2*M_PI) * delta_ktheta);
    for (len_t i = 0; i < num_ktheta; i++) {
        for (len_t j = 0; j < num_k; j++) {
            len_t idx = i * num_k + j;
            real_t dk_val = k_array[idx] - k0;
            real_t dkt_val = ktheta_array[idx] - ktheta0;
            amplitude[idx] = norm * std::exp(-0.5 * (dk_val*dk_val/(delta_k*delta_k) + 
                                                       dkt_val*dkt_val/(delta_ktheta*delta_ktheta)));
        }
    }
    
    std::cout << "WaveSpectrum: Set Gaussian spectrum centered at k=" << k0 
              << ", theta_k=" << ktheta0 << std::endl;
}

/**
 * Set custom spectrum from file
 */
void WaveSpectrum::setCustomSpectrum(const std::string& filename) {
    spectrum_type = OptionConstants::WAVE_SPECTRUM_CUSTOM;
    custom_file = filename;
    
    if (!loadCustomSpectrum()) {
        std::cerr << "WaveSpectrum: ERROR - Failed to load custom spectrum from " 
                  << filename << std::endl;
        return;
    }
    
    std::cout << "WaveSpectrum: Loaded custom spectrum from " << filename << std::endl;
}

/**
 * Load custom spectrum from file
 * Expected format: k, ktheta, amplitude (one mode per line, space or comma separated)
 */
bool WaveSpectrum::loadCustomSpectrum() {
    std::ifstream file(custom_file);
    if (!file.is_open()) {
        return false;
    }
    
    std::vector<real_t> k_temp, ktheta_temp, amp_temp;
    std::string line;
    
    while (std::getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Parse line (support both space and comma separation)
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream iss(line);
        
        real_t k_val, ktheta_val, amp_val;
        if (iss >> k_val >> ktheta_val >> amp_val) {
            k_temp.push_back(k_val);
            ktheta_temp.push_back(ktheta_val);
            amp_temp.push_back(amp_val);
        }
    }
    
    file.close();
    
    if (k_temp.size() != static_cast<size_t>(num_modes)) {
        std::cerr << "WaveSpectrum: WARNING - File contains " << k_temp.size() 
                  << " modes, but expected " << num_modes << std::endl;
        // Adjust grid sizes if needed
        // For now, we'll just use what we got
    }
    
    // Copy data to arrays
    num_modes = k_temp.size();
    num_k = std::sqrt(num_modes);  // Approximate
    num_ktheta = num_modes / num_k;
    
    k_array.assign(k_temp.begin(), k_temp.end());
    ktheta_array.assign(ktheta_temp.begin(), ktheta_temp.end());
    amplitude.assign(amp_temp.begin(), amp_temp.end());
    
    return true;
}

/**
 * Set amplitude for all modes (uniform amplitude)
 */
void WaveSpectrum::setAmplitude(real_t amp) {
    std::fill(amplitude.begin(), amplitude.end(), amp);
    
    std::cout << "WaveSpectrum: Set uniform amplitude = " << amp << std::endl;
}

/**
 * Set amplitude for a specific mode
 */
void WaveSpectrum::setAmplitude(len_t idx, real_t amp) {
    if (idx < num_modes) {
        amplitude[idx] = amp;
    } else {
        std::cerr << "WaveSpectrum: ERROR - Index " << idx 
                  << " out of range (max " << num_modes-1 << ")" << std::endl;
    }
}

/**
 * Set amplitudes from array
 */
void WaveSpectrum::setAmplitudes(const std::vector<real_t>& amps) {
    if (amps.size() != static_cast<size_t>(num_modes)) {
        std::cerr << "WaveSpectrum: WARNING - Amplitude array size (" 
                  << amps.size() << ") doesn't match num_modes (" 
                  << num_modes << ")" << std::endl;
        return;
    }
    
    amplitude.assign(amps.begin(), amps.end());
}

/**
 * Print spectrum information
 */
void WaveSpectrum::printInfo() const {
    std::cout << "=== Wave Spectrum Information ===" << std::endl;
    std::cout << "Spectrum type: ";
    switch (spectrum_type) {
        case OptionConstants::WAVE_SPECTRUM_UNIFORM:
            std::cout << "Uniform";
            break;
        case OptionConstants::WAVE_SPECTRUM_GAUSSIAN:
            std::cout << "Gaussian";
            break;
        case OptionConstants::WAVE_SPECTRUM_CUSTOM:
            std::cout << "Custom";
            break;
        default:
            std::cout << "Unknown";
    }
    std::cout << std::endl;
    
    std::cout << "Number of modes: " << num_modes << std::endl;
    std::cout << "Grid: " << num_k << " (k) x " << num_ktheta << " (theta_k)" << std::endl;
    
    if (!k_array.empty()) {
        auto minmax_k = std::minmax_element(k_array.begin(), k_array.end());
        auto minmax_kt = std::minmax_element(ktheta_array.begin(), ktheta_array.end());
        auto minmax_amp = std::minmax_element(amplitude.begin(), amplitude.end());
        
        std::cout << "k range: [" << *minmax_k.first << ", " << *minmax_k.second << "]" << std::endl;
        std::cout << "theta_k range: [" << *minmax_kt.first << ", " << *minmax_kt.second << "]" << std::endl;
        std::cout << "Amplitude range: [" << *minmax_amp.first << ", " << *minmax_amp.second << "]" << std::endl;
    }
    std::cout << "=================================" << std::endl;
}
