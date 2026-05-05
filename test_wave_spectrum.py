#!/usr/bin/env python3
"""
Test script for WaveSpectrum class.

This script tests the basic functionality of the WaveSpectrum class
before integrating it into the full DREAM simulation framework.
"""

import sys
sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM.Settings.WaveSpectrum import WaveSpectrum
import numpy as np


def test_uniform_spectrum():
    """Test uniform spectrum creation."""
    print("\n" + "="*60)
    print("TEST 1: Uniform Spectrum")
    print("="*60)
    
    ws = WaveSpectrum(num_k=50, num_ktheta=10)
    ws.setUniformSpectrum(
        k_min=35.0,
        k_max=45.0,
        ktheta_min=0.1,
        ktheta_max=0.3
    )
    ws.setAmplitude(1e-10)
    
    ws.printInfo()
    
    # Verify parameters
    assert ws.num_k == 50
    assert ws.num_ktheta == 10
    assert ws.num_modes == 500
    assert ws.spectrum_type == 'uniform'
    assert ws.amplitude == 1e-10
    
    print("✓ Uniform spectrum test PASSED\n")


def test_gaussian_spectrum():
    """Test Gaussian spectrum creation."""
    print("\n" + "="*60)
    print("TEST 2: Gaussian Spectrum")
    print("="*60)
    
    ws = WaveSpectrum(num_k=50, num_ktheta=10)
    ws.setGaussianSpectrum(
        k0=40.0,
        ktheta0=0.2,
        delta_k=2.0,
        delta_ktheta=0.05
    )
    ws.setAmplitude(5e-11)
    
    ws.printInfo()
    
    # Verify parameters
    assert ws.spectrum_type == 'gaussian'
    assert ws.parameters['k0'] == 40.0
    assert ws.parameters['delta_k'] == 2.0
    
    print("✓ Gaussian spectrum test PASSED\n")


def test_custom_spectrum():
    """Test custom spectrum setup (file loading not tested here)."""
    print("\n" + "="*60)
    print("TEST 3: Custom Spectrum (Setup Only)")
    print("="*60)
    
    ws = WaveSpectrum(num_k=100, num_ktheta=20)
    ws.setCustomSpectrum('/path/to/custom_spectrum.txt')
    ws.setAmplitude(1e-10)
    
    ws.printInfo()
    
    # Verify parameters
    assert ws.spectrum_type == 'custom'
    assert ws.parameters['filename'] == '/path/to/custom_spectrum.txt'
    
    print("✓ Custom spectrum setup test PASSED\n")


def test_individual_amplitudes():
    """Test setting individual amplitudes."""
    print("\n" + "="*60)
    print("TEST 4: Individual Amplitudes")
    print("="*60)
    
    ws = WaveSpectrum(num_k=10, num_ktheta=5)
    ws.setUniformSpectrum(35, 45, 0.1, 0.3)
    
    # Set different amplitudes for each mode
    amplitudes = np.linspace(1e-11, 1e-10, ws.num_modes)
    ws.setAmplitudes(amplitudes)
    
    ws.printInfo()
    
    # Verify
    assert ws.amplitudes_array is not None
    assert len(ws.amplitudes_array) == ws.num_modes
    assert np.isclose(np.min(ws.amplitudes_array), 1e-11)
    assert np.isclose(np.max(ws.amplitudes_array), 1e-10)
    
    print("✓ Individual amplitudes test PASSED\n")


def test_to_dict():
    """Test conversion to dictionary."""
    print("\n" + "="*60)
    print("TEST 5: Dictionary Conversion")
    print("="*60)
    
    ws = WaveSpectrum(num_k=20, num_ktheta=5)
    ws.setUniformSpectrum(35, 45, 0.1, 0.3)
    ws.setAmplitude(1e-10)
    
    d = ws.toDict()
    
    print("Dictionary keys:", list(d.keys()))
    print("num_k:", d['num_k'])
    print("num_ktheta:", d['num_ktheta'])
    print("spectrum_type:", d['spectrum_type'])
    print("amplitude:", d['amplitude'])
    
    # Verify
    assert 'num_k' in d
    assert 'num_ktheta' in d
    assert 'spectrum_type' in d
    assert 'parameters' in d
    assert 'amplitude' in d
    
    print("✓ Dictionary conversion test PASSED\n")


def main():
    """Run all tests."""
    print("\n" + "#"*60)
    print("# WaveSpectrum Class Test Suite")
    print("#"*60)
    
    try:
        test_uniform_spectrum()
        test_gaussian_spectrum()
        test_custom_spectrum()
        test_individual_amplitudes()
        test_to_dict()
        
        print("\n" + "#"*60)
        print("# ALL TESTS PASSED ✓")
        print("#"*60 + "\n")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
