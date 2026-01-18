"""Generates a LUKE equilibrium in HDF (.h5) format.

The magnetic field is constructed to mirror the parametrized "AnalyticB" field
used in DREAM unit testing.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).parents[2]
sys.path.append(str(ROOT / "tests" / "physics"))
sys.path.append(str(ROOT / "py"))

from numericmag import numericmag
from numericmag import savenummag

def main():
    # Generate the same analytic B that the Unittest
    # uses by default in UnitTest::InitializeAnalyticBRadialGridGenerator
    a = 2
    R0 = 4
    Delta_max = 0.6
    delta_max = 0.2
    g_min = 4
    g_max = 4.5
    kappa_min = 1.4
    kappa_max = 1.9
    delta_psi = np.pi

    nr = 50
    # introduce some safety margin in minor radius to avoid spline extrapolation
    r_profiles = np.linspace(0, 1.0*a, nr)

    def linear_profile(xmin, xmax):
        """Linear profile from xmin at rmin=0 to xmax at rmax=a.
        
        Linearly extrapolates for r beyond the minor radius a. 
        """
        rmin = 0
        rmax = a
        return xmin + (xmax - xmin) * (r_profiles - rmin) / (rmax - rmin)

    Delta = linear_profile(0, Delta_max)
    delta = linear_profile(0, delta_max)
    g_tor = linear_profile(g_min, g_max)
    kappa = linear_profile(kappa_min, kappa_max)
    psi_p0s = delta_psi * (r_profiles / a)**2 # reference poloidal flux profile
    psi_p0s = psi_p0s * R0  # use raw unnormalized psi

    nR = 40
    ntheta = 110
    mag = numericmag.constructMagneticField(
        Rp=R0, Zp=0, a=a, nR=nR, ntheta=ntheta,
        rG_R0=r_profiles, G_R0=g_tor, 
        rpsi=r_profiles, psi=psi_p0s,
        Delta=Delta, rDelta=r_profiles, 
        kappa=kappa, rkappa=r_profiles,
        delta=delta, rdelta=r_profiles, 
        retdict=True
    )

    filename = ROOT / "tests" / "data" / "analytic_b_luke_equilbrium.h5"
    savenummag.saveLUKE(str(filename), mag)

    r = np.linspace(0, a, nR+1)[1:]
    theta = np.linspace(0, 2*np.pi, ntheta+1)[:-1]
    mgR, mgT = np.meshgrid(r, theta)
    numericmag.plotMagneticField(mgR, mgT, mag["R"], mag["Z"], mag["Br"], mag["Bz"], mag["Bphi"])

if __name__ == "__main__":
    main()
