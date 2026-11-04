#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dream_pdrf_solver.py - A standalone, parameterized PDRF solver for DREAM.
Based on pdrf_SI.py by Hua-sheng XIE and Ziheng Zhou.
"""
import numpy as np
from scipy.linalg import eig as scipy_eig
import sympy as sp

class DreamPDRFSolver:
    """A parameterized PDRF solver for DREAM, based on pdrf_SI.py"""
    
    def __init__(self, B0=1.4, n_e=5e18, T_e=10.0, ion_mass_ratio=1836.0, Zeff=1.0):
        self.B0 = B0
        self.n_e = n_e
        self.T_e = T_e
        self.ion_mass_ratio = ion_mass_ratio
        self.Zeff = Zeff
        
        # Physical constants
        self.e_charge = 1.60217662e-19
        self.m_electron = 9.1094e-31
        self.epsilon0 = 8.854187817e-12
        self.c = 2.99792458e8
        
    def solve(self, kx, kz):
        """Solve dispersion relation for given kx, kz"""
        B0 = self.B0
        n_e = self.n_e
        T_e = self.T_e
        ion_mass_ratio = self.ion_mass_ratio
        Zeff = self.Zeff
        e_charge = self.e_charge
        m_electron = self.m_electron
        epsilon0 = self.epsilon0
        c = self.c
        
        # Species definition: [Electron, Ion]
        s = 2
        qs = np.array([-e_charge, Zeff * e_charge])
        ms = np.array([m_electron, ion_mass_ratio * m_electron])
        ns0 = np.array([n_e, n_e / Zeff])
        
        # Sound speeds (gamma = 5/3)
        gamma = 5.0 / 3.0
        T_joule = T_e * e_charge
        cs2 = np.array([
            gamma * T_joule / m_electron,
            gamma * T_joule / (ion_mass_ratio * m_electron)
        ])
        
        # Derived parameters
        wcs = qs * B0 / ms
        
        # Matrix dimensions
        dim = 4 * s + 6
        ind2 = 4 * s
        
        M = np.zeros((dim, dim), dtype=complex)
        A = np.eye(dim, dtype=complex)
        
        # Build matrices (Based on pdrf_SI.py logic)
        for j in range(s):
            ind = j * 4
            # dn~v
            M[0+ind, 1+ind] = -1j * kx * ns0[j]
            M[0+ind, 3+ind] = -1j * kz * ns0[j]
            
            # dv~n
            M[1+ind, 0+ind] = -1j * kx * cs2[j] / ns0[j]
            M[3+ind, 0+ind] = -1j * kz * cs2[j] / ns0[j]
            
            # dv~v (Cyclotron terms)
            M[1+ind, 2+ind] = wcs[j]
            M[2+ind, 1+ind] = -wcs[j]
            
            # dv~E
            M[1+ind, 0+ind2] = qs[j] / ms[j]
            M[2+ind, 1+ind2] = qs[j] / ms[j]
            M[3+ind, 2+ind2] = qs[j] / ms[j]
            
            # dE~v
            M[0+ind2, 1+ind] = -qs[j] * ns0[j] / epsilon0
            M[1+ind2, 2+ind] = -qs[j] * ns0[j] / epsilon0
            M[2+ind2, 3+ind] = -qs[j] * ns0[j] / epsilon0
            
        # EM Field terms
        M[0+ind2, 4+ind2] = -1j * kz * c**2
        M[1+ind2, 3+ind2] = 1j * kz * c**2
        M[1+ind2, 5+ind2] = -1j * kx * c**2
        M[2+ind2, 4+ind2] = 1j * kx * c**2
        
        M[3+ind2, 1+ind2] = 1j * kz
        M[4+ind2, 0+ind2] = -1j * kz
        M[4+ind2, 2+ind2] = 1j * kx
        M[5+ind2, 1+ind2] = -1j * kx
        
        try:
            # Use high-precision method matching pdrf_SI.py (oldmethod=2)
            # Step 1: Compute MA = A^{-1} * M
            MA = np.linalg.solve(A, M)
            
            # Step 2: Compute eigenvalues with standard numpy
            d_np = np.linalg.eigvals(MA)
            
            # Step 3: Apply high-precision conversion using sympy (16 digits)
            d = []
            for d_val in d_np:
                try:
                    d_sym = sp.N(d_val, 16)
                    d.append(complex(float(d_sym.real), float(d_sym.imag)))
                except Exception:
                    d.append(d_val)
            d = np.array(d, dtype=complex)
            
            # Convert eigenvalues: lambda = i*omega => omega = -i*lambda
            wtmp = 1j * d
            
            # Sort by real part (descending)
            inw1 = np.argsort(-np.real(wtmp))
            w1 = wtmp[inw1]
            
            return w1
        except Exception as e:
            print(f"DreamPDRFSolver Error: {e}")
            # Fallback to scipy method
            try:
                d, V = scipy_eig(M, A, check_finite=False)
                wtmp = 1j * d
                inw1 = np.argsort(-np.real(wtmp))
                return wtmp[inw1]
            except Exception as e2:
                print(f"Fallback method also failed: {e2}")
                return np.array([])

    def calculate_omega(self, k, theta):
        """Calculate omega for given k and theta"""
        kx = k * np.sin(theta)
        kz = k * np.cos(theta)
        return self.solve(kx, kz)
