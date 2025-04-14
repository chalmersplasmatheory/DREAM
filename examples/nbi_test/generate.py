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


#############################
# Simulation parameters     #
#############################

# Plasma parameters
T_initial = 500     # Initial electron temperature [eV]
Ip_initial = 400e5  # Initial plasma current [A]
TMAX = 1e-4         # Simulation time [s]

#Tokamak parameters
NR = 5              # Number of radial points
B0 = 5              # Magnetic field [T]
a = 0.5             # Plasma minor radius [m]
b = 0.6             # Wall radius [m]
R0 = 1.0            # Major radius [m]  

#NBI parameters
nbi_entry_point = "-1.5,-0.5,0.0" 
nbi_radius = 0.0775  
r_j_B = np.linspace(0, a, NR)  
j_B_values = 250e3 * np.ones(len(r_j_B))
nbi_direction_vector = "1.0,1.0,0.0"  # Direction of the NBI beam

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
#ds.radialgrid.setType(TYPE_ANALYTIC_TOROIDAL) 
#ds.radialgrid.setShaping(Delta=0.0, kappa=1.0, delta=0.0)  # Circular cross-section


# Set up time step
ds.timestep.setNt(10)
ds.timestep.setTmax(TMAX)
#ds.timestep.setIonization(dt0=4e-6, dtmax=4e-6, tmax=TMAX)
#ds.timestep.setMinSaveTimestep(3e-7)

# Add ion species
d = ds.eqsys.n_i.addIon('D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=1e19, T=500)
#ar = ds.eqsys.n_i.addIon('Ar', Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=2e18, T=500)


# Initial current / E-field setup
ds.eqsys.j_ohm.setInitialProfile(1, Ip0=Ip_initial)
ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0, R0=R0)

# Set up electron temperature and NBI
ds.eqsys.T_cold.setType(T_cold.TYPE_SELFCONSISTENT)
ds.eqsys.T_cold.setInitialProfile(T_initial)
ds.eqsys.T_cold.include_NBI = True


ds.eqsys.T_cold.setNBI({
   "ds": 0.1,
    "s_max": 4.05,
    "r_beam":  nbi_radius, 
    "P0": nbi_entry_point, 
    "n": nbi_direction_vector, 
    "Ti_beam": 4.8e-15,
    "m_i_beam": 3.344e-27,
    "beamPower": 1.6e6,
    "plasmaVolume": 2 * np.pi**2 * R0 * a**2,
    "Z0": 0,
    "Zion": 1,
    "R0": R0,
    "j_B/t": r_j_B,            
    "j_B/x": j_B_values,       
    "j_B/tinterp": 2    

})

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

def visualize_3d_tokamak():
    # Get NBI parameters
    P0 = np.array([float(x) for x in nbi_entry_point.split(',')])
    n = np.array([float(x) for x in "1.0,1.0,0.0".split(',')])
    s_max = 4.05
    r_beam = nbi_radius

    # Generate torus coordinates
    theta = np.linspace(0, 2*np.pi, 50)
    phi = np.linspace(0, 2*np.pi, 50)
    theta, phi = np.meshgrid(theta, phi)
    R = R0 + a * np.cos(theta)
    Z = a * np.sin(theta)
    X = R * np.cos(phi)
    Y = R * np.sin(phi)

    # Create 3D plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, color='royalblue', alpha=0.5, 
                   edgecolor='k', linewidth=0.1)

    # Add beam visualization
    s_vals = np.linspace(0, s_max, 30)
    phi_beam = np.linspace(0, 2*np.pi, 16)
    s_vals, phi_beam = np.meshgrid(s_vals, phi_beam)

    n_hat = n / np.linalg.norm(n)
    a_vec = np.array([0, 0, 1])
    if np.allclose(n_hat, a_vec):
        a_vec = np.array([0, 1, 0])
    e1 = a_vec - np.dot(a_vec, n_hat) * n_hat
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n_hat, e1)

    beam_center = P0[None, :] + s_vals[..., None] * n_hat
    beam_surface = (r_beam * np.cos(phi_beam)[..., None] * e1 +
                   r_beam * np.sin(phi_beam)[..., None] * e2 +
                   beam_center)

    Xb = beam_surface[..., 0]
    Yb = beam_surface[..., 1]
    Zb = beam_surface[..., 2]

    ax.plot_surface(Xb, Yb, Zb, color='orangered', alpha=0.7, 
                   edgecolor='k', linewidth=0.1)
    
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    ax.set_title('3D Tokamak and NBI Beam')
    ax.set_box_aspect([1,1,0.6])
    plt.show()

def visualize_flux_surfaces_top_view():
    # Create figure for top view
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Draw flux surfaces (circles in X-Y plane)
    radii = np.linspace(0, a, NR+1)
    theta = np.linspace(0, 2*np.pi, 100)
    
    # Create colormap
    max_index = NR-1
    colors = plt.cm.rainbow(np.linspace(0, 1, NR+2))
    
    ir = 0
    # Plot each flux surface
    for i, r in enumerate(radii):
        R1 = R0 + r * np.cos(0)
        R2 = R0 - r * np.cos(0)
        
        X1 = R1 * np.cos(theta)
        Y1 = R1 * np.sin(theta)
        X2 = R2 * np.cos(theta)
        Y2 = R2 * np.sin(theta)
        #ir = max_index - i
        # Use same color and label for matching inner/outer surfaces
        color_idx = ir
        label = f'ir = {ir-1}'
        ir= ir+1
        ax.plot(X1, Y1, '-', color=colors[ir], label=label)
        ax.plot(X2, Y2, '-', color=colors[ir])  # No label for second line
    
    # Add beam path
    P0 = np.array([float(x) for x in nbi_entry_point.split(',')])
    n = np.array([float(x) for x in "1.0,1.0,0.0".split(',')])
    s = np.linspace(0, 4.05, 100)
    
    beam_x = P0[0] + s * n[0]/np.sqrt(2)
    beam_y = P0[1] + s * n[1]/np.sqrt(2)
    ax.plot(beam_x, beam_y, 'r-', linewidth=2, label='NBI beam')
    
    # Add beam width
    n_perp = np.array([-n[1], n[0]]) / np.sqrt(2)
    for s_pos in np.linspace(0,4,200):
        pos_x = P0[0] + s_pos * n[0]/np.sqrt(2)
        pos_y = P0[1] + s_pos * n[1]/np.sqrt(2)
        width_x = pos_x + nbi_radius * np.array([-1, 1]) * n_perp[0]
        width_y = pos_y + nbi_radius * np.array([-1, 1]) * n_perp[1]
        ax.plot(width_x, width_y, 'r-', linewidth=2, alpha=0.5)

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title('Flux Surfaces and NBI Beam (Top View)')
    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

# Call both visualizations
visualize_3d_tokamak()
visualize_flux_surfaces_top_view()

print("Running DREAM simulation...")
ds.save('settings_nbi.h5')
runiface(ds, 'output_nbi.h5')


