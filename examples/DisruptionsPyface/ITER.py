import numpy as np
import sys

try:
    import DREAM
except ModuleNotFoundError:
    sys.path.append('~/Skrivbord/DREAM/py')
    import DREAM

from DREAM import DREAMSettings, runiface

import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Equations.IonSpecies as Ions

# ITER parameters
R0 = 6.                            # major radius (m)
a  = 2.                            # minor radius (m)
b  = 2.833                         # minor radius of tokamak wall (m)
B0 = 5.3                           # toroidal magnetic field on-axis (T)

# Timescales
tmax_TQ_max = 2e-2                 # Maximum time for temporary thermal quench simulation (s)
tmax_CQ_max = 4e-1                 # Total simulation time (s)

# Simulation parameters
ne0     = 1e20                     # electron density (m^-3)
T_core  = 2e4                      # Temperature at plasma core (eV)
T_edge  = 2e2                      # Temperature at plasma edge (eV)
tau_w   = 0.5                      # Wall time (s)
Ip      = 15e6                     # Target plasma current (A)
Ire_tol = 150e3                    # Maximum RE current tolerated (A)
j_exp   = 0.41                     # Plasma current profile parameter for exponent
j_scale = (1 - 0.001**(1./j_exp))  # Plasma current profile parameter for scaling

# Radial grid resolution
NR = 101


integrator = None

def setRadialGrid(ds, nr=NR, toroidal=True):
    """
    Set an ITER-like magnetic field for the radial grid generator.

    Source: István Pusztai.
    """
    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)

    if toroidal:
        ds.radialgrid.setMajorRadius(R0)
    else:
        ds.radialgrid.setMajorRadius(np.inf)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(b)
    ds.radialgrid.setNr(nr)

    # Radial grid for analytical magnetic field
    r = np.linspace(0, a, nr)

    def f(a0, a1, a2, a3, a4):
        return np.polyval(np.flip([a0, a1, a2, a3, a4]), r)

    # Shaping parameters
    ds.radialgrid.setShapeParameter('kappa', r=r, data=f(1.5, 0, 0, 0, .02))
    ds.radialgrid.setShapeParameter('delta', r=r, data=f(0, 0, .035, 0, .017))
    ds.radialgrid.setShapeParameter('Delta', r=r, data=f(0, 0, -.00658678, -.00634124, 0))

    # Poloidal magnetic flux
    if toroidal:
        ds.radialgrid.setShapeParameter('psi_p0', r=r, data=f(0, 0, .794102, -0.117139, 0))
    else:
        ds.radialgrid.setShapeParameter('psi_p0', r=r, data=np.zeros(r.shape))

    # Toroidal magnetic field function
    ds.radialgrid.setShapeParameter('GOverR0', r=r, data=f(B0, 0, -.310741, .107719, 0))

def getInitialTemperatureProfile(T_core, T_edge, nr=NR):
    """  Returns the initial temperature profile. """
    r = np.linspace(0, a, nr)
    T = T_core * (1 - (1 - T_edge/T_core) * (r/a)**2)
    return r, T

def getInitialCurrentDensityProfile(j_scale, j_exp, Ip, nr=NR):
    """ Returns the initial current density profile. """
    r = np.linspace(0, a, nr)
    j = (1 - j_scale*(r/a)**2)**j_exp
    return r, j, Ip

def getIntegrator(nr=NR):
    global integrator
    if integrator is None:
        ds_integrate = DREAMSettings()
        ds_integrate.hottailgrid.setEnabled(False)
        ds_integrate.runawaygrid.setEnabled(False)
        setRadialGrid(ds_integrate, nr=nr)
        ds_integrate.eqsys.n_i.addIon('D', n=ne0, Z=1, Z0=1, iontype=Ions.IONS_DYNAMIC, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
        ds_integrate.eqsys.E_field.setPrescribedData(0.001)
        ds_integrate.eqsys.T_cold.setPrescribedData(temperature=20)
        ds_integrate.timestep.setNt(1)
        ds_integrate.timestep.setTmax(1)
        integrator = runiface(ds_integrate, quiet=True)
    return integrator, integrator.grid.r