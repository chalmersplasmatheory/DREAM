"""Generates a LUKE equilibrium in HDF (.h5) format.

The magnetic field is constructed to mirror the parametrized "AnalyticB" field
used in DREAM unit testing.
"""
from pathlib import Path
import sys

import numpy as np
from scipy.constants import mu_0

ROOT = Path(__file__).parents[2]
sys.path.append(str(ROOT / "tests" / "physics"))
sys.path.append(str(ROOT / "py"))

from numericmag import numericmag
from numericmag import savenummag


def generate(
    index, a, R0, G, kappa, psi, delta, Delta,
    nR=40, ntheta=110, nRref=50, plot=False
):
    """
    Generate a magnetic field.
    """
    rref = np.linspace(0, 1.0*a, nRref)
    rho = rref / a

    nR = 40
    ntheta = 110
    mag = numericmag.constructMagneticField(
        Rp=R0, Zp=0, a=a, nR=nR, ntheta=ntheta,
        rG_R0=rref, G_R0=G(rho), 
        rpsi=rref, psi=psi(rho),
        Delta=Delta(rho), rDelta=rref, 
        kappa=kappa(rho), rkappa=rref,
        delta=delta(rho), rdelta=rref, 
        retdict=True
    )

    filename = ROOT / "tests" / "data" / f"analytic_b_luke_equilbrium_{index:d}.h5"
    savenummag.saveLUKE(str(filename), mag)

    if plot:
        r = np.linspace(0, a, nR+1)[1:]
        theta = np.linspace(0, 2*np.pi, ntheta+1)[:-1]
        mgR, mgT = np.meshgrid(r, theta)
        numericmag.plotMagneticField(mgR, mgT, mag["R"], mag["Z"], mag["Br"], mag["Bz"], mag["Bphi"])


def main():
    plot = False

    def linear_profile(rho, xmin, xmax):
        return xmin + (xmax-xmin) * rho

    # Ola's equilibrium
    generate(
        index=0, a=2, R0=4,
        G=lambda r : linear_profile(r, 4, 4.5),
        kappa=lambda r : linear_profile(r, 1.4, 1.9),
        psi=lambda r : 4 * np.pi * r**2,
        Delta=lambda r : linear_profile(r, 0, 0.6),
        delta=lambda r : linear_profile(r, 0, 0.2),
        plot=plot
    )

    # Mathias' multiple Bmin equilibrium
    IpRef = 11e6
    generate(
        index=1, a=1, R0=4,
        G=lambda r : (4.0 + 0.1*r),
        kappa=lambda r : 1+1.9*r**3,
        psi=lambda r : -mu_0 * IpRef * (1-r**3) * 1,
        Delta=lambda r : 0.05*r,
        delta=lambda r : -0.3*np.sqrt(r),
        plot=plot
    )


if __name__ == "__main__":
    sys.exit(main())


