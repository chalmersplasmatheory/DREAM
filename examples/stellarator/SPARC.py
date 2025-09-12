import numpy as np
import sys
import netCDF4 
sys.path.append('$DREAMPATH/py')
from DREAM import DREAMIO
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d as interp1d

import DREAM

import DREAM.Settings.RadialGrid as RadialGrid


shaping = DREAMIO.LoadHDF5AsDict('SPARCdata/shaping.h5')
transp = netCDF4.Dataset('SPARCdata/10000.CDF','r'); # TRANSP

time_slice 	= 5.973
time 		= transp.variables['TIME'][:]
idx_time 	= np.argmin(np.abs(time-time_slice))

# SPARC parameters
# Major radius
R0_t = transp.variables['RAXIS'][idx_time] / 100             # major radius (m)
R = transp.variables['RMAJM'][idx_time,:] / 100
idx_R 		= np.argmin(np.abs(R-R0_t))
r_data = (R[idx_R+1:] - R0_t)

a  = r_data[-1]          # minor radius (m)
b  = 1.183*a            # minor radius of tokamak wall (m)
B0 = shaping['GOverR0'][0]           # toroidal magnetic field on-axis (T)
Ip = float(transp.variables['PCUR'][idx_time])          # Target plasma current (A)
R0 = shaping['R0'][()]

# Timescales
tmax_TQ_max = 2e-3     # Maximum time for temporary thermal quench simulation (s)
tmax_CQ = 20e-3         # Total simulation time (s)

# Simulation parameters
ne_data = transp.variables['NE'][idx_time,:]*100**3  # electron density (m^-3)
Te_data = transp.variables['TE'][idx_time,:]
T_initial = Te_data[0]

johm_data = transp.variables['CUR'][idx_time,:] * 100**2
johm_data /= np.max(johm_data)

ne0 = ne_data[0]

tau_w = 80e-3

# Radial grid resolution
NR = 101

def setRadialGrid(ds, nr=NR):
    """
    Set an SPARC-like magnetic field for the radial grid generator.
    """
    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)

    ds.radialgrid.setMajorRadius(R0)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(b)
    ds.radialgrid.setNr(nr)

    r = np.linspace(0, a, (nr-1)*2+1)[1:]

    r_d = shaping['r']

    psi = interp1d(r_d, shaping['psi'], kind='cubic', fill_value="extrapolate")(r)
    kappa = interp1d(r_d,shaping['kappa'], kind='linear', fill_value="extrapolate")(r)
    delta = interp1d(r_d, savgol_filter(shaping['delta'], 10, 2), kind='cubic', fill_value="extrapolate")(r)
    Delta = interp1d(r_d,shaping['Delta'], kind='cubic', fill_value="extrapolate")(r)
    Delta[0] = Delta[1]
    GOverR0 = interp1d(r_d,shaping['GOverR0'], kind='cubic', fill_value="extrapolate")(r)

    # Shaping parameters 
    ds.radialgrid.setShapeParameter('kappa', r=r, data=kappa)
    ds.radialgrid.setShapeParameter('delta', r=r, data=delta)
    ds.radialgrid.setShapeParameter('Delta', r=r, data=Delta)

    # Poloidal magnetic flux
    ds.radialgrid.setShapeParameter('psi_p0', r=r, data=psi)

    # Toroidal magnetic field function
    ds.radialgrid.setShapeParameter('GOverR0', r=r, data=GOverR0)

def getInitialTemperatureProfile(nr=NR):
    """  Returns the initial temperature profile. """
    r = np.linspace(0, a, (nr-1)*20+1)[1:]
    T = np.interp(r, r_data, Te_data)
    return r, T

def getInitialDensityProfile(nr=NR):
    """  Returns the initial temperature profile. """
    r = np.linspace(0, a, (nr-1)*20+1)[1:]
    n = np.interp(r, r_data, ne_data)
    return r, n

def getInitialCurrentDensityProfile(nr=NR):
    """ Returns the initial current density profile. """
    r = np.linspace(0, a, (nr-1)*20+1)[1:]
    j = np.interp(r, r_data, johm_data)
    return r, j