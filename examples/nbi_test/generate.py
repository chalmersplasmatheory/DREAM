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

#Tokamak parameters - ROME
NR = 50              # Number of radial points
B0 = 5              # Magnetic field [T]
a = 0.23             # Plasma minor radius [m]
b = 0.23           # Wall radius [m]
R0 =  0.798            # Major radius [m]  


#NBI parameters TO COMPARE WITH ROME
r= np.linspace(0,a,NR+1)  # Radial grid points
n_i_profile = 3e19 * (1 - (r/a)**2)
nbi_entry_point = [0.685,-1.028,0.0] 
nbi_radius = 0.0775  
nbi_direction_vector = [0.0,1.0,0.0] 


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


Ti_profile = 250   # constant 250 eV
# Add ionized deuterium
ds.eqsys.n_i.addIon(name='D_core', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=n_i_profile, r=r, T=Ti_profile)
#ds.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=v, T=1)


# Initial current / E-field setup
ds.eqsys.j_ohm.setInitialProfile(1, Ip0=Ip_initial)
ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0, R0=R0)


nbi = NBISettings()
nbi.setEnabled(True)
##ROME
#nbi.setTCVGaussian(True)  # Use TCV Gaussian beam profile for TCV
nbi.setOrigin(nbi_entry_point)
nbi.setDirection(nbi_direction_vector)

nbi.setBeamParameters(r_beam=nbi_radius, Ti_beam=25*1.6021e-16, m_i_beam=3.344e-27)
nbi.setIons(Z0=0, Zion=1)  
nbi.setPower(beam_power=1) #57e3
nbi.setR0_NBI(R0)

nbi.visualize_flux_surfaces_top_view(r)



##TCV SETTINGS
#nbi.setPower(beam_power=1.03e6)
#nbi.setBeamParameters(Ti_beam= 28*1.6021e-16)
#nbi.setTCVGaussian(True)  # Use TCV Gaussian beam profile for TCV
#nbi.visualize_3d_tokamak(a = 0.23)
#nbi.visualize_flux_surfaces_top_view(np.linspace(0,a,20))


#t_beam = np.linspace(0,tMax, 100)
#P_beam = np.zeros_like(t_beam)
#P_beam[(t_beam >= 0) & (t_beam <= float(NBI_duriation) - 1e-10)] = beam_power
#nbi.setPowerProfile(t_beam, P_beam)





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


