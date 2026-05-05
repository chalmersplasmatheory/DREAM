#ifndef _DREAM_SETTINGS_WAVE_SPECTRUM_HPP
#define _DREAM_SETTINGS_WAVE_SPECTRUM_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include <vector>
#include <string>

namespace DREAM {
    /**
     * Class for managing wave spectrum parameters for quasilinear diffusion.
     * 
     * This class stores the discretized wave spectrum in (k, theta_k) space
     * and provides methods to set different types of spectra (uniform, Gaussian, custom).
     */
    class WaveSpectrum {
    private:
        len_t num_k;           // Number of wavenumber grid points
        len_t num_ktheta;      // Number of angle grid points
        len_t num_modes;       // Total number of modes = num_k * num_ktheta
        
        std::vector<real_t> k_array;       // Wavenumber array [num_modes]
        std::vector<real_t> ktheta_array;  // Angle array [num_modes]
        std::vector<real_t> amplitude;     // Wave amplitude for each mode [num_modes]
        
        enum OptionConstants::wave_spectrum_type spectrum_type;
        
        // Parameters for different spectrum types
        struct UniformParams {
            real_t k_min, k_max;
            real_t ktheta_min, ktheta_max;
        } uniform_params;
        
        struct GaussianParams {
            real_t k0, ktheta0;
            real_t delta_k, delta_ktheta;
        } gaussian_params;
        
        std::string custom_file;  // File path for custom spectrum
        
    public:
        /**
         * Constructor
         * 
         * @param nk       Number of wavenumber grid points
         * @param nkt      Number of angle grid points
         */
        WaveSpectrum(len_t nk=100, len_t nkt=20);
        
        /**
         * Destructor
         */
        ~WaveSpectrum();
        
        /**
         * Set uniform spectrum in (k, theta_k) space
         * 
         * @param k_min        Minimum wavenumber
         * @param k_max        Maximum wavenumber
         * @param ktheta_min   Minimum angle
         * @param ktheta_max   Maximum angle
         */
        void setUniformSpectrum(real_t k_min, real_t k_max, 
                               real_t ktheta_min, real_t ktheta_max);
        
        /**
         * Set Gaussian spectrum centered at (k0, ktheta0)
         * 
         * @param k0           Center wavenumber
         * @param ktheta0      Center angle
         * @param delta_k      Wavenumber spread
         * @param delta_ktheta Angle spread
         */
        void setGaussianSpectrum(real_t k0, real_t ktheta0, 
                                real_t delta_k, real_t delta_ktheta);
        
        /**
         * Set custom spectrum from file
         * 
         * @param filename  Path to file containing spectrum data
         *                  Format: k, ktheta, amplitude (one mode per line)
         */
        void setCustomSpectrum(const std::string& filename);
        
        /**
         * Set amplitude for all modes (uniform amplitude)
         * 
         * @param amp  Wave amplitude (normalized to n_e m_e c^2)
         */
        void setAmplitude(real_t amp);
        
        /**
         * Set amplitude for a specific mode
         * 
         * @param idx  Mode index
         * @param amp  Wave amplitude
         */
        void setAmplitude(len_t idx, real_t amp);
        
        /**
         * Set amplitudes from array
         * 
         * @param amps  Array of amplitudes [num_modes]
         */
        void setAmplitudes(const std::vector<real_t>& amps);
        
        // Accessor methods
        len_t getNumK() const { return num_k; }
        len_t getNumKtheta() const { return num_ktheta; }
        len_t getNumModes() const { return num_modes; }
        
        real_t getK(len_t idx) const { return k_array[idx]; }
        real_t getKtheta(len_t idx) const { return ktheta_array[idx]; }
        real_t getAmplitude(len_t idx) const { return amplitude[idx]; }
        
        const std::vector<real_t>& getKArray() const { return k_array; }
        const std::vector<real_t>& getKthetaArray() const { return ktheta_array; }
        const std::vector<real_t>& getAmplitudeArray() const { return amplitude; }
        
        enum OptionConstants::wave_spectrum_type getSpectrumType() const { 
            return spectrum_type; 
        }
        
        /**
         * Check if spectrum is properly initialized
         */
        bool isInitialized() const { return !k_array.empty(); }
        
        /**
         * Print spectrum information
         */
        void printInfo() const;
        
    private:
        /**
         * Initialize arrays with given sizes
         */
        void initializeArrays();
        
        /**
         * Load custom spectrum from file
         */
        bool loadCustomSpectrum();
    };
}

#endif /* _DREAM_SETTINGS_WAVE_SPECTRUM_HPP */
