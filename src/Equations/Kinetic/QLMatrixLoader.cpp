#include "DREAM/Equations/Kinetic/QLMatrixLoader.hpp"
#include <iostream>
#include <stdexcept>
#include <cstring>

#ifdef BUILD_WITH_HDF5

namespace DREAM {
namespace Equations {
namespace Kinetic {

QLMatrixLoader::QLMatrixLoader(const std::string& hdf5_file)
    : filename(hdf5_file),
      k_array(nullptr), ktheta_array(nullptr), omega_array(nullptr),
      num_modes(0),
      p_array(nullptr), xi_array(nullptr),
      num_p(0), num_xi(0),
      D_pp_base(nullptr), D_pxi_base(nullptr), D_xixi_base(nullptr),
      A0(1.0), A_current(1.0)
{
    loadFromHDF5();
}

QLMatrixLoader::~QLMatrixLoader() {
    if (k_array) delete[] k_array;
    if (ktheta_array) delete[] ktheta_array;
    if (omega_array) delete[] omega_array;
    if (p_array) delete[] p_array;
    if (xi_array) delete[] xi_array;
    if (D_pp_base) delete[] D_pp_base;
    if (D_pxi_base) delete[] D_pxi_base;
    if (D_xixi_base) delete[] D_xixi_base;
}

void QLMatrixLoader::loadFromHDF5() {
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        throw std::runtime_error("Failed to open HDF5 file: " + filename);
    }
    
    std::cerr << "Loading quasi-linear diffusion matrix from " << filename << "..." << std::endl;
    
    // Read metadata
    hid_t attr_id = H5Aopen(file_id, "A0", H5P_DEFAULT);
    if (attr_id >= 0) {
        H5Aread(attr_id, H5T_NATIVE_DOUBLE, &A0);
        H5Aclose(attr_id);
    } else {
        std::cerr << "Warning: Could not read A0 attribute, using default value 1.0" << std::endl;
        A0 = 1.0;
    }
    
    // Read spectrum grid
    {
        hid_t dataset_id = H5Dopen(file_id, "spectrum/k", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: spectrum/k");
        }
        
        hid_t space_id = H5Dget_space(dataset_id);
        hsize_t dims[1];
        H5Sget_simple_extent_dims(space_id, dims, nullptr);
        num_modes = static_cast<len_t>(dims[0]);
        
        k_array = new real_t[num_modes];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, k_array);
        H5Dclose(dataset_id);
        H5Sclose(space_id);
    }
    
    {
        hid_t dataset_id = H5Dopen(file_id, "spectrum/ktheta", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: spectrum/ktheta");
        }
        
        ktheta_array = new real_t[num_modes];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ktheta_array);
        H5Dclose(dataset_id);
    }
    
    {
        hid_t dataset_id = H5Dopen(file_id, "spectrum/omega", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: spectrum/omega");
        }
        
        omega_array = new real_t[num_modes];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, omega_array);
        H5Dclose(dataset_id);
    }
    
    // Read momentum grid
    {
        hid_t dataset_id = H5Dopen(file_id, "momentum/p", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: momentum/p");
        }
        
        hid_t space_id = H5Dget_space(dataset_id);
        hsize_t dims[1];
        H5Sget_simple_extent_dims(space_id, dims, nullptr);
        num_p = static_cast<len_t>(dims[0]);
        
        p_array = new real_t[num_p];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_array);
        H5Dclose(dataset_id);
        H5Sclose(space_id);
    }
    
    {
        hid_t dataset_id = H5Dopen(file_id, "momentum/xi", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: momentum/xi");
        }
        
        hid_t space_id = H5Dget_space(dataset_id);
        hsize_t dims[1];
        H5Sget_simple_extent_dims(space_id, dims, nullptr);
        num_xi = static_cast<len_t>(dims[0]);
        
        xi_array = new real_t[num_xi];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xi_array);
        H5Dclose(dataset_id);
        H5Sclose(space_id);
    }
    
    // Read base diffusion tensors
    {
        hid_t dataset_id = H5Dopen(file_id, "diffusion/D_pp_base", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: diffusion/D_pp_base");
        }
        
        D_pp_base = new real_t[num_p * num_xi];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D_pp_base);
        H5Dclose(dataset_id);
    }
    
    {
        hid_t dataset_id = H5Dopen(file_id, "diffusion/D_pxi_base", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: diffusion/D_pxi_base");
        }
        
        D_pxi_base = new real_t[num_p * num_xi];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D_pxi_base);
        H5Dclose(dataset_id);
    }
    
    {
        hid_t dataset_id = H5Dopen(file_id, "diffusion/D_xixi_base", H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: diffusion/D_xixi_base");
        }
        
        D_xixi_base = new real_t[num_p * num_xi];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D_xixi_base);
        H5Dclose(dataset_id);
    }
    
    H5Fclose(file_id);
    
    std::cerr << "✓ Loaded:" << std::endl;
    std::cerr << "  Wave modes: " << num_modes << std::endl;
    std::cerr << "  Momentum grid: " << num_p << " × " << num_xi << std::endl;
    std::cerr << "  Reference amplitude A0: " << A0 << std::endl;
    
    // Count non-zero entries for verification
    len_t nnz_pp = 0, nnz_pxi = 0, nnz_xixi = 0;
    for (len_t idx = 0; idx < num_p * num_xi; idx++) {
        if (D_pp_base[idx] != 0.0) nnz_pp++;
        if (D_pxi_base[idx] != 0.0) nnz_pxi++;
        if (D_xixi_base[idx] != 0.0) nnz_xixi++;
    }
    std::cerr << "  Non-zero D_pp entries: " << nnz_pp << std::endl;
    std::cerr << "  Non-zero D_pξ entries: " << nnz_pxi << std::endl;
    std::cerr << "  Non-zero D_ξξ entries: " << nnz_xixi << std::endl;
}

real_t QLMatrixLoader::getDppBase(len_t i, len_t j) const {
    if (i >= num_p || j >= num_xi) {
        std::cerr << "Error: Index out of bounds in getDppBase (" 
                  << i << ", " << j << ")" << std::endl;
        return 0.0;
    }
    return D_pp_base[i * num_xi + j];
}

real_t QLMatrixLoader::getDpxiBase(len_t i, len_t j) const {
    if (i >= num_p || j >= num_xi) {
        std::cerr << "Error: Index out of bounds in getDpxiBase (" 
                  << i << ", " << j << ")" << std::endl;
        return 0.0;
    }
    return D_pxi_base[i * num_xi + j];
}

real_t QLMatrixLoader::getDxixiBase(len_t i, len_t j) const {
    if (i >= num_p || j >= num_xi) {
        std::cerr << "Error: Index out of bounds in getDxixiBase (" 
                  << i << ", " << j << ")" << std::endl;
        return 0.0;
    }
    return D_xixi_base[i * num_xi + j];
}

}}} // namespace DREAM::Equations::Kinetic

#else // !BUILD_WITH_HDF5

// Stub implementation when HDF5 is not available
namespace DREAM {
namespace Equations {
namespace Kinetic {

QLMatrixLoader::QLMatrixLoader(const std::string&) {
    throw std::runtime_error("HDF5 support not enabled in DREAM build. "
                             "Please rebuild with HDF5 support to use pre-computed matrices.");
}

QLMatrixLoader::~QLMatrixLoader() {}

void QLMatrixLoader::loadFromHDF5() {
    throw std::runtime_error("HDF5 support not enabled");
}

real_t QLMatrixLoader::getDppBase(len_t, len_t) const { return 0; }
real_t QLMatrixLoader::getDpxiBase(len_t, len_t) const { return 0; }
real_t QLMatrixLoader::getDxixiBase(len_t, len_t) const { return 0; }

}}} // namespace DREAM::Equations::Kinetic

#endif // BUILD_WITH_HDF5
