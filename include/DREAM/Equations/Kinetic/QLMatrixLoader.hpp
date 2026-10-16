#ifndef _DREAM_EQUATIONS_KINETIC_QL_MATRIX_LOADER_HPP
#define _DREAM_EQUATIONS_KINETIC_QL_MATRIX_LOADER_HPP

#include "DREAM/config.h"
#include <string>

#ifdef BUILD_WITH_HDF5
#include <hdf5.h>
#endif

namespace DREAM {
namespace Equations {
namespace Kinetic {

/**
 * Loads pre-computed quasi-linear diffusion matrix from HDF5 file.
 * 
 * The matrix structure is time-independent and computed by Python.
 * C++ applies time-dependent amplitude scaling during runtime.
 * 
 * The diffusion tensor is factorized as:
 *   D_ij(p, ξ, t) = D_ij_base(p, ξ) * |A(t)|² / |A₀|²
 * 
 * where D_ij_base is loaded from HDF5 (computed with reference amplitude A₀),
 * and A(t) is the time-dependent wave amplitude set at runtime.
 */
class QLMatrixLoader {
private:
    std::string filename;
    
    // Loaded data
    real_t *k_array;           ///< Wavenumber array [num_modes]
    real_t *ktheta_array;      ///< Wave angle array [num_modes]
    real_t *omega_array;       ///< Frequency array [num_modes]
    len_t num_modes;           ///< Number of wave modes
    
    real_t *p_array;           ///< Momentum array [num_p]
    real_t *xi_array;          ///< Pitch-angle cosine array [num_xi]
    len_t num_p;               ///< Number of momentum grid points
    len_t num_xi;              ///< Number of pitch-angle grid points
    
    real_t *D_pp_base;         ///< Base D_pp tensor [num_p * num_xi]
    real_t *D_pxi_base;        ///< Base D_pξ tensor [num_p * num_xi]
    real_t *D_xixi_base;       ///< Base D_ξξ tensor [num_p * num_xi]
    
    // Reference amplitude used in Python pre-computation
    real_t A0;
    
    // Current time-dependent amplitude (set by C++ at runtime)
    real_t A_current;

public:
    /**
     * Constructor: loads data from HDF5 file
     * @param hdf5_file  Path to the HDF5 file containing pre-computed matrix
     */
    QLMatrixLoader(const std::string& hdf5_file);
    
    /**
     * Destructor: frees allocated memory
     */
    ~QLMatrixLoader();
    
    /**
     * Load data from HDF5 file
     * @throws std::runtime_error if file cannot be opened or read
     */
    void loadFromHDF5();
    
    // ========================================================================
    // Getters for grid information
    // ========================================================================
    
    len_t getNumModes() const { return num_modes; }
    len_t getNumP() const { return num_p; }
    len_t getNumXi() const { return num_xi; }
    
    real_t getK(len_t m) const { return k_array[m]; }
    real_t getKtheta(len_t m) const { return ktheta_array[m]; }
    real_t getOmega(len_t m) const { return omega_array[m]; }
    
    real_t getP(len_t i) const { return p_array[i]; }
    real_t getXi(len_t j) const { return xi_array[j]; }
    
    // ========================================================================
    // Access to base diffusion coefficients (per unit amplitude²)
    // ========================================================================
    
    /**
     * Get base D_pp coefficient at grid point (i, j)
     * @param i  Momentum index
     * @param j  Pitch-angle index
     * @return D_pp_base[i, j]
     */
    real_t getDppBase(len_t i, len_t j) const;
    
    /**
     * Get base D_pξ coefficient at grid point (i, j)
     * @param i  Momentum index
     * @param j  Pitch-angle index
     * @return D_pξ_base[i, j]
     */
    real_t getDpxiBase(len_t i, len_t j) const;
    
    /**
     * Get base D_ξξ coefficient at grid point (i, j)
     * @param i  Momentum index
     * @param j  Pitch-angle index
     * @return D_ξξ_base[i, j]
     */
    real_t getDxixiBase(len_t i, len_t j) const;
    
    // ========================================================================
    // Time-dependent amplitude control
    // ========================================================================
    
    /**
     * Set current wave amplitude (time-dependent)
     * @param A_t  Current amplitude at time t
     */
    void setCurrentAmplitude(real_t A_t) { A_current = A_t; }
    
    /**
     * Get current wave amplitude
     * @return Current amplitude A(t)
     */
    real_t getCurrentAmplitude() const { return A_current; }
    
    /**
     * Get reference amplitude used in Python pre-computation
     * @return Reference amplitude A₀
     */
    real_t getReferenceAmplitude() const { return A0; }
    
    // ========================================================================
    // Scaled diffusion coefficients (with time-dependent amplitude)
    // ========================================================================
    
    /**
     * Get scaled D_pp coefficient at time t
     * @param i  Momentum index
     * @param j  Pitch-angle index
     * @return D_pp(i, j, t) = D_pp_base[i, j] * |A(t)|² / |A₀|²
     */
    real_t getDppScaled(len_t i, len_t j) const {
        return D_pp_base[i * num_xi + j] * (A_current * A_current) / (A0 * A0);
    }
    
    /**
     * Get scaled D_pξ coefficient at time t
     * @param i  Momentum index
     * @param j  Pitch-angle index
     * @return D_pξ(i, j, t) = D_pξ_base[i, j] * |A(t)|² / |A₀|²
     */
    real_t getDpxiScaled(len_t i, len_t j) const {
        return D_pxi_base[i * num_xi + j] * (A_current * A_current) / (A0 * A0);
    }
    
    /**
     * Get scaled D_ξξ coefficient at time t
     * @param i  Momentum index
     * @param j  Pitch-angle index
     * @return D_ξξ(i, j, t) = D_ξξ_base[i, j] * |A(t)|² / |A₀|²
     */
    real_t getDxixiScaled(len_t i, len_t j) const {
        return D_xixi_base[i * num_xi + j] * (A_current * A_current) / (A0 * A0);
    }
};

}}} // namespace DREAM::Equations::Kinetic

#endif /* _DREAM_EQUATIONS_KINETIC_QL_MATRIX_LOADER_HPP */
