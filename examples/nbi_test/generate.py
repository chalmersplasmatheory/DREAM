#!/usr/bin/env python3
import sys
sys.path.append('../../py')


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from DREAM import DREAMSettings, runiface
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Solver as Solver
from DREAM.Settings.Equations import IonSpecies
import DREAM.Settings.Equations.IonSpecies as Ions
from DREAM.Settings.Equations.NBISettings import NBISettings 


#############################
# Simulation parameters     #
#############################

# Plasma parameters
T_initial = 800     # Initial electron temperature [eV]
Ip_initial = 120e5  # Initial plasma current [A]
TMAX = 1e-6         # Simulation time [s]

#Tokamak parameters
NR = 20              # Number of radial points
B0 = 5              # Magnetic field [T]
a = 0.23             # Plasma minor radius [m]
b = 0.23           # Wall radius [m]
R0 =  0.798            # Major radius [m]  

#NBI parameters
nbi_entry_point = [0.685,-1.028,0.0]
#nbi_entry_point = [0.585,-1.028,0.0]
#nbi_entry_point = [0,-1,0.0]
nbi_radius = 0.07  # Beam radius [m]  
r_j_B = np.linspace(0, nbi_radius, 10)  
j_B_values = 250e3 * np.ones(len(r_j_B))
#nbi_direction_vector = [1.0,1.0,0.0]  # Direction of the NBI beam
nbi_direction_vector = [0.0,1.0,0.0]  # Direction of the NBI beam
s_max = 2

#############################
# Initialize DREAM settings #
#############################
ds = DREAMSettings()

# Set up grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(b)
ds.radialgrid.setNr(NR) 
ds.radialgrid.setMajorRadius(R0)


# Set up time step
ds.timestep.setNt(1)
ds.timestep.setTmax(TMAX)
r = np.linspace(0, a, NR*10)
n_profile = 3e19*(1-(r/(a))**2) 


# Add ion species
# Fully ionized ions (core plasma)
n_i_profile = 3e19 * (1 - (r/a)**2)  # peak = 3e19 m⁻³
Ti_profile = 250   # constant 250 eV

# Neutrals (background gas)
n_neutral_profile = 7e14 * np.ones_like(r) # peak = 7e14 m⁻³ try parabilic
Tn_profile = 250        # same temp for simplicity obs
# Add ionized deuterium
ds.eqsys.n_i.addIon(name='D_core', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=n_i_profile, r=r, T=Ti_profile)

# Add neutral deuterium
#ds.eqsys.n_i.addIon(name='D_neutral', Z=1, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_neutral_profile, r=r, T=Tn_profile)


# Initial current / E-field setup
ds.eqsys.j_ohm.setInitialProfile(1, Ip0=Ip_initial)
ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0, R0=R0)

nbi = NBISettings()
nbi.setEnabled(True)
nbi.setOrigin(nbi_entry_point)
nbi.setDirection(nbi_direction_vector)
nbi.setBeamParameters(r_beam=nbi_radius, Ti_beam=4.8e-15, m_i_beam=3.344e-27, s_max=s_max)
nbi.setIons(Z0=0, Zion=1)  # Fully ionized deuterium
nbi.setPower(beam_power=1.0)
nbi.setCurrentProfile(j_B_t=r_j_B, j_B_x=j_B_values, tinterp=0)
nbi.setR0_NBI(R0)
nbi.visualize_3d_tokamak(nbi_entry_point, nbi_direction_vector, nbi_radius, R0, a, s_max)
nbi.visualize_flux_surfaces_top_view(nbi_entry_point, nbi_direction_vector, nbi_radius, R0, a, s_max, NR)


# Set up electron temperature and NBI
ds.eqsys.T_cold.setType(T_cold.TYPE_SELFCONSISTENT)
ds.eqsys.T_cold.setInitialProfile(T_initial)
ds.eqsys.T_cold.setNBI(nbi)

print(ds.eqsys.T_cold.todict())
print("include_NBI setting is:", ds.eqsys.T_cold.include_NBI)


# Runaway electron settings
ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

# Disable kinetic grids
ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

# Solver settings
ds.solver.setType(Solver.LINEAR_IMPLICIT)
ds.solver.setVerbose(True)
# Include fluid
ds.other.include('fluid')







print("Running DREAM simulation...")
ds.save('settings_nbi.h5')
runiface(ds, 'output_nbi.h5')


