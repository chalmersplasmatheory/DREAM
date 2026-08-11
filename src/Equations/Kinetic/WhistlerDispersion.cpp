/**
 * WhistlerDispersion class - Simplified version for reference.
 * 
 * NOTE: This class is primarily kept for:
 * 1. Simplified dispersion relation calculations (used in Python pre-computation)
 * 2. Documentation of physical parameters and constants
 * 3. Historical reference for PDRF integration approach
 * 
 * Current workflow uses Python pre-computation with BON solver,
 * not runtime PDRF calculation. See:
 * - /data/zhzhou/DREAM/examples/test_dispersion_relation/precompute_ql_matrix.py
 * - QLMatrixLoader for C++ HDF5 loading
 */

// ============================================================================
// DEPRECATED: Python C API + PDRF Integration
// ============================================================================
// The following code is kept for historical reference but is NOT USED
// in the current pre-computed matrix workflow.
// Current approach: Python pre-computation with BON solver → HDF5 → C++ load
// ============================================================================

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "DREAM/Equations/Kinetic/WhistlerDispersion.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <sstream>

using namespace DREAM;

// Global Python objects (initialized once)
// WARNING: These are NOT USED in current pre-computed workflow
static PyObject* pdrf_module = nullptr;
static PyObject* pdrf_solver_class = nullptr;
static bool python_initialized = false;

/**
 * Initialize Python interpreter and import PDRF module
 * DEPRECATED: Not used in pre-computed matrix workflow
 */
static bool initializePython() {
    if (python_initialized) {
        return true;
    }
    
    // Initialize Python interpreter
    Py_Initialize();
    
    if (!Py_IsInitialized()) {
        std::cerr << "Error: Failed to initialize Python interpreter" << std::endl;
        return false;
    }
    
    // Import NumPy C API
    import_array1(false);  // Required for NumPy C API
    
    // Add DREAM PDRF path to sys.path
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.insert(0, '/data/zhzhou/DREAM/include/DREAM/Equations/Kinetic')");
    
    // Import DREAM PDRF module (integrated version)
    std::cerr << "[PDRF] Attempting to import dream_pdrf..." << std::endl;
    pdrf_module = PyImport_ImportModule("dream_pdrf");
    if (!pdrf_module) {
        std::cerr << "[PDRF] ERROR: Failed to import dream_pdrf module" << std::endl;
        PyErr_Print();  // Print Python error traceback
        Py_Finalize();
        return false;
    }
    std::cerr << "[PDRF] ✓ Successfully imported dream_pdrf" << std::endl;
    
    // Get DreamPDRFSolver class
    pdrf_solver_class = PyObject_GetAttrString(pdrf_module, "DreamPDRFSolver");
    if (!pdrf_solver_class || !PyCallable_Check(pdrf_solver_class)) {
        std::cerr << "Error: Cannot find DreamPDRFSolver class" << std::endl;
        Py_DECREF(pdrf_module);
        Py_Finalize();
        return false;
    }
    
    python_initialized = true;
    std::cout << "✓ Python PDRF backend initialized successfully" << std::endl;
    
    return true;
}

/**
 * Cleanup Python resources
 */
static void cleanupPython() {
    if (python_initialized) {
        Py_XDECREF(pdrf_solver_class);
        Py_XDECREF(pdrf_module);
        Py_Finalize();
        python_initialized = false;
    }
}

/**
 * Constructor
 * 
 * NOTE: In current pre-computed workflow, this class is primarily used for:
 * - Calculating physical constants (ω_ce, ω_pe, v_A, w_factor)
 * - Simplified dispersion relation (when use_simple_dispersion=true)
 * 
 * The PDRF integration code below is NOT USED in production.
 */
WhistlerDispersion::WhistlerDispersion(real_t B0_val, real_t density_val, 
                                      real_t ion_mass_ratio_val, 
                                      real_t Zeff_val,
                                      bool use_simple)
    : B0(B0_val), density(density_val), 
      ion_mass_ratio(ion_mass_ratio_val), Zeff(Zeff_val),
      use_simple_dispersion(use_simple) {
    
    initializeParameters();
    
    // Initialize Python on first instance (only if not using simple method)
    if (!use_simple_dispersion && !initializePython()) {
        std::cerr << "Warning: Python initialization failed, dispersion calculations will not work" << std::endl;
    }
}

/**
 * Destructor
 */
WhistlerDispersion::~WhistlerDispersion() {
    // Note: We don't call cleanupPython() here because other instances might exist
    // Python will be cleaned up when the program exits
}

/**
 * Initialize derived parameters
 */
void WhistlerDispersion::initializeParameters() {
    // Calculate cyclotron frequencies: ω_c = qB/m
    omega_ce = e_charge * B0 / m_electron;  // Electron cyclotron frequency
    omega_ci = omega_ce / ion_mass_ratio;   // Ion cyclotron frequency
    
    // Calculate plasma frequencies: ω_p = sqrt(n q^2 / (ε₀ m))
    omega_pe = std::sqrt(density * e_charge * e_charge / (epsilon0 * m_electron));
    omega_pi = omega_pe / std::sqrt(ion_mass_ratio);
    
    // Calculate Alfvén velocity: v_A = B / sqrt(μ₀ n_i m_i)
    real_t m_ion = ion_mass_ratio * m_electron;
    real_t n_ion = density / Zeff;  // Quasi-neutrality: n_e = Zeff * n_i
    v_A = B0 / std::sqrt(mu0 * n_ion * m_ion);
    
    // Calculate simplified dispersion coefficient: w = ω_ce * c² / ω_pe²
    // This is the correct cold plasma whistler approximation from BON code
    w_factor = omega_ce * c * c / (omega_pe * omega_pe);
    
    if (use_simple_dispersion) {
        std::cout << "WhistlerDispersion initialized (SIMPLIFIED dispersion relation):" << std::endl;
        std::cout << "  B0 = " << B0 << " T" << std::endl;
        std::cout << "  n_e = " << density << " m^-3" << std::endl;
        std::cout << "  v_A = " << v_A << " m/s" << std::endl;
        std::cout << "  w = ω_ce * c² / ω_pe² = " << w_factor << " m²/s" << std::endl;
        std::cout << "  Dispersion: ω = k|k_∥| * w" << std::endl;
    } else {
        std::cout << "WhistlerDispersion initialized (Python C API + PDRF):" << std::endl;
        std::cout << "  B0 = " << B0 << " T" << std::endl;
        std::cout << "  n_e = " << density << " m^-3" << std::endl;
        std::cout << "  ω_ce = " << omega_ce << " rad/s (" 
                  << omega_ce/(2*M_PI)/1e9 << " GHz)" << std::endl;
        std::cout << "  ω_pe = " << omega_pe << " rad/s (" 
                  << omega_pe/(2*M_PI)/1e9 << " GHz)" << std::endl;
        std::cout << "  ω_ci = " << omega_ci << " rad/s (" 
                  << omega_ci/(2*M_PI)/1e6 << " MHz)" << std::endl;
    }
}

/**
 * Call Python PDRF solver to calculate wave frequency
 * 
 * DEPRECATED: This method is NOT USED in the current pre-computed matrix workflow.
 * Current approach: Python pre-computation with BON solver → HDF5 → C++ load via QLMatrixLoader
 * 
 * This code is kept for historical reference only.
 */
std::vector<real_t> WhistlerDispersion::solvePDRFMethod(real_t kx, real_t kz) const {
    if (!python_initialized) {
        std::cerr << "Error: Python not initialized" << std::endl;
        return {};
    }
    
    // Create DreamPDRFSolver instance with parameters
    // Constructor: DreamPDRFSolver(B0, n_e, T_e=10.0)
    PyObject* solver_args = PyTuple_Pack(3,
                                         PyFloat_FromDouble(B0),
                                         PyFloat_FromDouble(density),
                                         PyFloat_FromDouble(10.0));  // T_e = 10 eV
    
    PyObject* solver_instance = PyObject_CallObject(pdrf_solver_class, solver_args);
    if (!solver_instance) {
        PyErr_Print();
        std::cerr << "Error: Failed to create DreamPDRFSolver instance" << std::endl;
        Py_DECREF(solver_args);
        return {};
    }
    Py_DECREF(solver_args);
    
    // For now, use a simpler approach: call solve method with kx, kz
    PyObject* solve_method = PyObject_GetAttrString(solver_instance, "solve");
    if (!solve_method || !PyCallable_Check(solve_method)) {
        std::cerr << "Error: Cannot find solve method" << std::endl;
        Py_DECREF(solver_instance);
        return {};
    }
    
    // Call solve(kx, kz) - simplified interface
    PyObject* args = PyTuple_Pack(2, 
                                  PyFloat_FromDouble(kx),
                                  PyFloat_FromDouble(kz));
    
    PyObject* result = PyObject_CallObject(solve_method, args);
    
    if (!result) {
        PyErr_Print();
        std::cerr << "Error: PDRF solve() failed" << std::endl;
        Py_DECREF(args);
        Py_DECREF(solve_method);
        Py_DECREF(solver_instance);
        return {};
    }
    
    // Parse result: dream_pdrf returns numpy array of eigenvalues directly
    if (!PyArray_Check(result)) {
        std::cerr << "Error: Invalid result format (expected numpy array)" << std::endl;
        Py_DECREF(result);
        Py_DECREF(args);
        Py_DECREF(solve_method);
        Py_DECREF(solver_instance);
        return {};
    }
    
    // Extract eigenvalues from numpy array
    PyArrayObject* w1_array = reinterpret_cast<PyArrayObject*>(result);
    npy_intp n_freqs = PyArray_SIZE(w1_array);
    
    std::vector<real_t> frequencies;
    frequencies.reserve(n_freqs);
    
    for (npy_intp i = 0; i < n_freqs; i++) {
        // Get complex value
        char* ptr = static_cast<char*>(PyArray_DATA(w1_array)) + i * PyArray_ITEMSIZE(w1_array);
        
        // For complex128, extract real and imaginary parts
        double re, im;
        memcpy(&re, ptr, sizeof(double));
        memcpy(&im, ptr + sizeof(double), sizeof(double));
        
        // Check if valid: positive real part, small imaginary part
        if (re > 0 && std::abs(im) / std::abs(re) < 1e-6) {
            frequencies.push_back(re);
        }
    }
    
    // Cleanup
    Py_DECREF(result);
    Py_DECREF(args);
    Py_DECREF(solve_method);
    Py_DECREF(solver_instance);
    
    // std::cout << "  Found " << frequencies.size() << " valid frequencies via Python C API" << std::endl;
    
    return frequencies;
}

/**
 * Calculate wave frequency omega for given (k, theta)
 * 
 * CURRENT USAGE:
 * - If use_simple_dispersion=true: Uses simplified formula ω = k|k_∥| * w
 *   (This IS used by Python pre-computation script)
 * 
 * DEPRECATED PATH:
 * - If use_simple_dispersion=false: Would call PDRF solver (NOT USED)
 * 
 * For full dispersion calculations, use Python pre-computation with BON solver instead.
 */
real_t WhistlerDispersion::calculateOmega(real_t k, real_t theta) const {
    // If using simplified dispersion, call the simple method
    if (use_simple_dispersion) {
        return calculateOmegaSimple(k, theta);
    }
    
    // Otherwise use full PDRF calculation
    // Calculate wave vector components
    // Following PDRF convention: B0 along z-axis, k in x-z plane
    real_t kx = k * std::sin(theta);  // Perpendicular component
    real_t kz = k * std::cos(theta);  // Parallel component
    
    // std::cout << "\nCalculating dispersion for k=" << k << " m^-1, θ=" 
    //           << theta*180/M_PI << "°" << std::endl;
    // std::cout << "  k_x = " << kx << " m^-1, k_z = " << kz << " m^-1" << std::endl;
    
    // Solve via Python PDRF (C API - fast!)
    std::vector<real_t> frequencies = solvePDRFMethod(kx, kz);
    
    if (frequencies.empty()) {
        std::cerr << "Error: No valid frequencies found from PDRF solver" << std::endl;
        return -1.0;
    }
    
    // Select whistler branch: Ω_ci << ω < Ω_ce
    real_t omega_selected = -1.0;
    for (real_t omega : frequencies) {
        if (omega > 10*omega_ci && omega < 0.5*omega_ce) {
            omega_selected = omega;
            break;  // Take first valid whistler frequency
        }
    }
    
    if (omega_selected < 0) {
        std::cerr << "Warning: No whistler branch found, taking lowest frequency" << std::endl;
        std::sort(frequencies.begin(), frequencies.end());
        omega_selected = frequencies[0];
    }
    
    // std::cout << "  Selected ω = " << omega_selected << " rad/s ("
    //           << omega_selected/(2*M_PI)/1e6 << " MHz)" << std::endl;
    // std::cout << "  ω/ω_ci = " << omega_selected/omega_ci << std::endl;
    // std::cout << "  ω/ω_ce = " << omega_selected/omega_ce << std::endl;
    
    return omega_selected;
}

/**
 * Calculate group velocity dω/dk
 */
real_t WhistlerDispersion::calculateGroupVelocity(real_t k, real_t theta) const {
    const real_t dk = 1e-3 * k;  // Small perturbation
    
    real_t omega1 = calculateOmega(k - dk/2, theta);
    real_t omega2 = calculateOmega(k + dk/2, theta);
    
    if (omega1 < 0 || omega2 < 0) {
        std::cerr << "Error: Failed to calculate group velocity" << std::endl;
        return 0.0;
    }
    
    return (omega2 - omega1) / dk;
}

/**
 * Calculate polarization vectors
 */
void WhistlerDispersion::calculatePolarization(real_t omega, real_t k, real_t theta,
                                              real_t &Ex_out, real_t &Ey_out, real_t &Ez_out) const {
    // TODO: Implement polarization calculation from eigenvectors
    Ex_out = 1.0;
    Ey_out = 0.0;
    Ez_out = 0.0;
}

/**
 * Calculate normalization factor
 */
real_t WhistlerDispersion::calculateDenominator(real_t omega, real_t k, real_t theta) const {
    // TODO: Implement proper denominator calculation
    return 1.0;
}

/**
 * Check if frequency is in whistler range
 */
bool WhistlerDispersion::isInWhistlerRange(real_t omega) const {
    return (omega > 10*omega_ci && omega < 0.5*omega_ce);
}

/**
 * Print dispersion relation information
 */
void WhistlerDispersion::printInfo() const {
    std::cout << "\n=== Whistler Dispersion Relation Info ===" << std::endl;
    std::cout << "B0 = " << B0 << " T" << std::endl;
    std::cout << "n_e = " << density << " m^-3" << std::endl;
    std::cout << "ω_ce = " << omega_ce/(2*M_PI)/1e9 << " GHz" << std::endl;
    std::cout << "ω_pe = " << omega_pe/(2*M_PI)/1e9 << " GHz" << std::endl;
    std::cout << "ω_ci = " << omega_ci/(2*M_PI)/1e6 << " MHz" << std::endl;
    if (use_simple_dispersion) {
        std::cout << "Method: Simplified whistler dispersion (ω = k|k_∥| * w)" << std::endl;
    } else {
        std::cout << "Method: Python C API + PDRF (embedded)" << std::endl;
    }
    std::cout << "=====================================\n" << std::endl;
}

/**
 * Calculate wave frequency using simplified whistler dispersion relation.
 * 
 * For Ω_ci ≪ ω ≪ Ω_ce, the cold plasma whistler dispersion simplifies to:
 *   ω = k|k_∥| * w,  where w = ω_ce * c² / ω_pe²
 * 
 * This is derived from the cold plasma dielectric tensor in the whistler limit.
 * Reference: calculate_kperp_bon.py (BON code)
 */
real_t WhistlerDispersion::calculateOmegaSimple(real_t k, real_t theta) const {
    // Calculate parallel wavenumber: k_∥ = k * cos(θ)
    real_t k_parallel = k * std::cos(theta);
    
    // Simplified whistler dispersion: ω = k|k_∥| * w
    real_t omega = k * std::abs(k_parallel) * w_factor;
    
    return omega;
}

/**
 * Calculate perpendicular wavenumber using simplified whistler dispersion.
 * 
 * Inverse problem: given ω and k_∥, find k_⊥ using:
 *   ω = k|k_∥| * w  =>  k = ω / (|k_∥| * w)
 *   k_⊥ = sqrt(k² - k_∥²)
 * 
 * Reference: calculate_kperp_bon.py calculate_kperp_resonance()
 */
real_t WhistlerDispersion::calculateKperpSimple(real_t omega, real_t k_par) const {
    // Avoid division by zero
    if (std::abs(k_par) < 1e-12) {
        std::cerr << "Warning: k_∥ is too small for simplified dispersion calculation" << std::endl;
        return -1.0;
    }
    
    // Calculate total wavenumber: k = ω / (|k_∥| * w)
    real_t k_total = omega / (std::abs(k_par) * w_factor);
    
    // Calculate k_⊥ = sqrt(k² - k_∥²)
    real_t k_perp_sq = k_total * k_total - k_par * k_par;
    
    if (k_perp_sq < 0) {
        // No real solution - wave cannot propagate at this frequency
        return -1.0;
    }
    
    return std::sqrt(k_perp_sq);
}
