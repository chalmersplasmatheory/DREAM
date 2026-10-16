"""
WaveSpectrum class for DREAM quasilinear diffusion settings.

This module provides a Python interface for configuring wave spectra
used in quasilinear diffusion calculations.
"""

import numpy as np


class WaveSpectrum:
    """
    Class for managing wave spectrum parameters in quasilinear diffusion.
    
    This class stores the discretized wave spectrum in (k, theta_k) space
    and provides methods to set different types of spectra.
    """
    
    def __init__(self, num_k=100, num_ktheta=20):
        """
        Constructor.
        
        Args:
            num_k (int): Number of wavenumber grid points
            num_ktheta (int): Number of angle grid points
        """
        self.num_k = num_k
        self.num_ktheta = num_ktheta
        self.num_modes = num_k * num_ktheta
        
        self.spectrum_type = None
        self.parameters = {}
        self.amplitude = None
        self.amplitudes_array = None
        
    def setUniformSpectrum(self, k_min, k_max, ktheta_min, ktheta_max):
        """
        Set uniform spectrum in (k, theta_k) space.
        
        Args:
            k_min (float): Minimum wavenumber
            k_max (float): Maximum wavenumber
            ktheta_min (float): Minimum angle
            ktheta_max (float): Maximum angle
        """
        self.spectrum_type = 'uniform'
        self.parameters = {
            'k_min': k_min,
            'k_max': k_max,
            'ktheta_min': ktheta_min,
            'ktheta_max': ktheta_max
        }
        
        print(f"WaveSpectrum: Set uniform spectrum with "
              f"{self.num_k} x {self.num_ktheta} = {self.num_modes} modes")
        print(f"  k range: [{k_min}, {k_max}]")
        print(f"  theta_k range: [{ktheta_min}, {ktheta_max}]")
        
    def setGaussianSpectrum(self, k0, ktheta0, delta_k, delta_ktheta):
        """
        Set Gaussian spectrum centered at (k0, ktheta0).
        
        Args:
            k0 (float): Center wavenumber
            ktheta0 (float): Center angle
            delta_k (float): Wavenumber spread (standard deviation)
            delta_ktheta (float): Angle spread (standard deviation)
        """
        self.spectrum_type = 'gaussian'
        self.parameters = {
            'k0': k0,
            'ktheta0': ktheta0,
            'delta_k': delta_k,
            'delta_ktheta': delta_ktheta
        }
        
        # Generate grid covering ±3 sigma
        k_min = max(k0 - 3 * delta_k, 1e-10)
        k_max = k0 + 3 * delta_k
        ktheta_min = max(ktheta0 - 3 * delta_ktheta, 1e-10)
        ktheta_max = ktheta0 + 3 * delta_ktheta
        
        # Store the actual grid ranges
        self.parameters['k_min'] = k_min
        self.parameters['k_max'] = k_max
        self.parameters['ktheta_min'] = ktheta_min
        self.parameters['ktheta_max'] = ktheta_max
        
        print(f"WaveSpectrum: Set Gaussian spectrum centered at "
              f"k={k0}, theta_k={ktheta0}")
        
    def setCustomSpectrum(self, filename):
        """
        Set custom spectrum from file.
        
        Args:
            filename (str): Path to file containing spectrum data
                           Format: k, ktheta, amplitude (one mode per line)
        """
        self.spectrum_type = 'custom'
        self.parameters = {'filename': filename}
        
        print(f"WaveSpectrum: Will load custom spectrum from {filename}")
        
    def setAmplitude(self, amplitude):
        """
        Set amplitude for all modes (uniform amplitude).
        
        Args:
            amplitude (float): Wave amplitude (normalized to n_e m_e c^2)
        """
        self.amplitude = amplitude
        self.amplitudes_array = None
        
        print(f"WaveSpectrum: Set uniform amplitude = {amplitude}")
        
    def setAmplitudes(self, amplitudes):
        """
        Set amplitudes for each mode individually.
        
        Args:
            amplitudes (array-like): Array of amplitudes [num_modes]
        """
        amplitudes = np.asarray(amplitudes)
        if len(amplitudes) != self.num_modes:
            raise ValueError(f"Amplitude array size ({len(amplitudes)}) "
                           f"doesn't match num_modes ({self.num_modes})")
        
        self.amplitudes_array = amplitudes
        self.amplitude = None
        
        print(f"WaveSpectrum: Set individual amplitudes for {self.num_modes} modes")
        
    def toDict(self):
        """
        Convert to dictionary for HDF5 output.
        
        Returns:
            dict: Dictionary representation of wave spectrum
        """
        result = {
            'num_k': self.num_k,
            'num_ktheta': self.num_ktheta,
            'spectrum_type': self.spectrum_type,
            'parameters': self.parameters
        }
        
        if self.amplitude is not None:
            result['amplitude'] = self.amplitude
        elif self.amplitudes_array is not None:
            result['amplitudes'] = self.amplitudes_array.tolist()
        
        return result
    
    def printInfo(self):
        """Print spectrum information."""
        print("=== Wave Spectrum Information ===")
        print(f"Spectrum type: {self.spectrum_type}")
        print(f"Number of modes: {self.num_modes}")
        print(f"Grid: {self.num_k} (k) x {self.num_ktheta} (theta_k)")
        
        if self.parameters:
            print(f"Parameters: {self.parameters}")
        
        if self.amplitude is not None:
            print(f"Amplitude: {self.amplitude}")
        elif self.amplitudes_array is not None:
            print(f"Amplitude range: [{np.min(self.amplitudes_array)}, "
                  f"{np.max(self.amplitudes_array)}]")
        
        print("=================================")
