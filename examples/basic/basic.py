#!/usr/bin/env python3
#
# A very basic DREAM Python example. This script generates a basic
# DREAM input file which can be passed to 'dreami'.
#
# Run as
#
#   $ ./basic.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')
from DREAMSettings import DREAMSettings

import Equations.IonSpecies as Ions
import Solver


ds = DREAMSettings()

times  = [0]
radius = [0, 1]

# Set E_field
efield = 1e-7*np.ones((len(times), len(radius)))
ds.equationsystem.E_field.setPrescribedData(efield=efield, times=times, radius=radius)

# Set n_cold
density = 1e19 * np.ones((len(times), len(radius)))
ds.equationsystem.n_cold.setPrescribedData(density=density, times=times, radius=radius)

# Set temperature
temperature = 1000 * np.ones((len(times), len(radius)))
ds.equationsystem.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)

# Set ions
ds.equationsystem.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=5e19)
ds.equationsystem.n_i.addIon(name='He', Z=2, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e19)
ds.equationsystem.n_i.addIon(name='B', Z=5, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=.2e19)
ds.equationsystem.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=8e19)



# Hot-tail grid settings
pmax = 5
ds.hottailgrid.setNxi(30)
ds.hottailgrid.setNp(100)
ds.hottailgrid.setPmax(pmax)

# Set initial hot electron distribution function
fhot_r = np.array([0])
fhot_p = np.linspace(0, pmax, 100)
fhot_xi = np.array([1])
nR, nP, nXi = fhot_r.size, fhot_p.size, fhot_xi.size
fhot = np.zeros((nR, nXi, nP))
for k in range(0, nR):
    for j in range(0, nXi):
        fhot[k,j,:] = (pmax - fhot_p) / pmax

ds.equationsystem.f_hot.setInitialValue(init=fhot, r=fhot_r, p=fhot_p, xi=fhot_xi)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(10)

# Use the linear solver
ds.solver.setType(Solver.TYPE_LINEAR_IMPLICIT)

# Set time stepper
ds.timestep.setTmax(1.0)
ds.timestep.setNt(4)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

