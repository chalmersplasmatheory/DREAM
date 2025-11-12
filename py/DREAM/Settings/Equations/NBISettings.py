import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants for NBI Gaussian profile types
NBI_PROFILE_TCV = 1
NBI_PROFILE_ITER = 2
NBI_PROFILE_CUSTOM = 3
    
class NBISettings:
    """
    Handles Neutral Beam Injection settings in the DREAM framework.
    This class is used by ColdElectronTemperature to manage NBI-specific parameters.
    Pre-set values match TCV parameters.
    """

    
    def __init__(self):
        self.enabled = False
        self.s_max = 2.05 # Maximum beam path length [m]
        self.r_beam = 0.089 # Beam radius [m]
        self.P0 = [-1.2,-0.42,0.0] # Beam entry point [m]
        self.n = [0.5547, 0.83205, 0.0] # Beam direction vector (unit vector)
        self.Ti_beam = 28*1.6021e-16 # Beam energy [J]
        self.m_i_beam = 3.344e-27 # Beam ion mass [kg] (Deuterium)
        self.R0 = 0.88 # Major radius [m]
        self.j_B_t = np.linspace(0, 0.0775, 50) # Beam current profile time points [s] (Matching ROME radius)
        self.j_B_x = 250e3 * np.ones(len(self.j_B_t)) # Beam current profile values [A/m^2]
        self.j_B_tinterp = 0 # Interpolation method for time profile
        self.gaussian_profile = 0 # 0=disabled, 1=TCV, 2=ITER
        self.a = 0.23 # Plasma minor radius [m]
        self.P_NBI_t = None # Beam power profile time points [s]
        self.P_NBI_x = None # Beam power profile values [W]
        self.P_NBI_tinterp =0 # Interpolation method for power profile
        self.energy_fractions = [1,0,0]# Fractions of beam energy for multi-energy components

    def setEnabled(self, enabled=True):
        """Enable/disable NBI."""
        self.enabled = enabled  
        
    def setCurrentProfile(self, profile_type, j_B_t=None, j_B_x=None, tinterp=0):
        """Set NBI Gaussian profile type. profile_type: 0=disabled, 1=TCV, 2=ITER, 3 = custom"""
        if profile_type not in [1,2,3,0]:
            raise ValueError("Invalid profile_type. 1 (TCV), 2 (ITER), or 3 (custom).")
        
        elif  profile_type == 3:
            if j_B_t is None or j_B_x is None:
                raise ValueError("For custom profile, j_B_t and j_B_x must be provided.")
            self.j_B_t = j_B_t
            self.j_B_x = j_B_x
            self.j_B_tinterp = tinterp

        self.gaussian_profile = profile_type
        
    def setOrigin(self, P0):
        """Set beam entry point."""
        self.P0 = P0

    def setDirection(self, n):
        """Set beam direction vector."""
        self.n = n
        
    def setEnergyFractions(self, f_full=None, f_half=None, f_third=None):
        """Set energy fractions for multi-energy beam components. NBI beams contain a mix of atomic (full energy) and molecular species that
        dissociate into half and third energy components. """
        fractions = [f_full, f_half, f_third]
        if abs(sum(fractions) - 1.0) > 1e-6:
            raise ValueError(f"Energy fractions must sum to 1.0, got {sum(fractions):.4f}")

        self.energy_fractions = fractions

    def setBeamParameters(self, r_beam=None, Ti_beam=None, m_i_beam=None, s_max=None):
        """Set beam physical parameters."""
        if r_beam is not None: self.r_beam = r_beam
        if Ti_beam is not None: self.Ti_beam = Ti_beam
        if m_i_beam is not None: self.m_i_beam = m_i_beam
        if s_max is not None: self.s_max = s_max

    def setRadius(self, R0=None, a=None):
        if R0 is not None: self.R0 = R0
        if a is not None: self.a = a

    def setPowerProfile(self, P_NBI_t = None, P_NBI_x = None, tinterp=0):
        """Set beam power profile in one dimension."""
        if P_NBI_x is not None:
            if np.isscalar(P_NBI_x):
                # Convert scalar to array with single value
                self.P_NBI_x = [float(P_NBI_x)]
                if P_NBI_t is None:
                    self.P_NBI_t = [0.0] 
                else: 
                    raise ValueError("Warning: P_NBI_t ignored when P_NBI_x is scalar.")
            else:
                self.P_NBI_x = P_NBI_x
                if P_NBI_t is not None:
                    self.P_NBI_t = P_NBI_t 
                else:
                    raise ValueError("Warning: P_NBI_t must be provided when P_NBI_x is an array.")
    
        else: 
            raise ValueError("NBI: P_NBI_x must be provided.")
            
        self.P_NBI_tinterp = tinterp

    def visualize_3d_tokamak(self):
        """
        Visualize the tokamak as a 3D torus and the NBI beam as a cylindrical surface.
        """
        P0 = np.array(self.P0)
        n = np.array(self.n)
        # Create torus surface
        theta, phi = np.meshgrid(np.linspace(0, 2*np.pi, 50), np.linspace(0, 2*np.pi, 50))
        R = self.R0 + self.a * np.cos(theta)
        Z = self.a * np.sin(theta)
        X = R * np.cos(phi)
        Y = R * np.sin(phi)

        # Create 3D figure
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, color='royalblue', alpha=0.5, edgecolor='k', linewidth=0.1)

        # Beam unit vector and orthogonal basis
        n_hat = n / np.linalg.norm(n)
        a_vec = np.array([0, 0, 1]) if not np.allclose(n_hat, [0, 0, 1]) else np.array([0, 1, 0])
        e1 = a_vec - np.dot(a_vec, n_hat) * n_hat
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n_hat, e1)

        # Generate beam cylinder surface
        s_vals, phi_beam = np.meshgrid(np.linspace(0, self.s_max, 30), np.linspace(0, 2*np.pi, 16))
        beam_center = P0[None, :] + s_vals[..., None] * n_hat
        beam_surface = (self.r_beam * np.cos(phi_beam)[..., None] * e1 +
                        self.r_beam * np.sin(phi_beam)[..., None] * e2 +
                        beam_center)

        ax.plot_surface(beam_surface[..., 0], beam_surface[..., 1], beam_surface[..., 2],
                        color='orangered', alpha=0.7, edgecolor='k', linewidth=0.1)

        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
        ax.set_title('3D Tokamak and NBI Beam')
        ax.set_box_aspect([1, 1, 0.6])
        plt.tight_layout()
        plt.show()


    def visualize_flux_surfaces_top_view(self, radius_vector):
        """
        Visualize circular flux surfaces and NBI beam in the top-down (X-Y) view.
        """
        fig, ax = plt.subplots(figsize=(10, 10))
        NR =len(radius_vector)
        r_f = radius_vector
        r_mid = 0.5 * (r_f[:-1] + r_f[1:])
        theta = np.linspace(0, 2*np.pi, 100)
        colors = plt.cm.rainbow(np.linspace(0, 1, NR + 2))

        for i, r in enumerate(r_mid):
            R1 = self.R0 + r
            R2 = self.R0 - r
            X1 = R1 * np.cos(theta)
            Y1 = R1 * np.sin(theta)
            X2 = R2 * np.cos(theta)
            Y2 = R2 * np.sin(theta)
            ax.plot(X1, Y1, '-', color=colors[i], label=f'ir = {i}')
            ax.plot(X2, Y2, '-', color=colors[i])

        # Beam path
        P0 = np.array(self.P0)
        n = np.array(self.n)
        n= n / np.linalg.norm(n)
        s = np.linspace(0, self.s_max, 100)
        beam_x = P0[0] + s * n[0]
        beam_y = P0[1] + s * n[1]
        ax.plot(beam_x, beam_y, 'r-', linewidth=2, label='NBI beam')

        # Beam width lines
        n_perp = np.array([-n[1], n[0]]) / np.linalg.norm([n[0], n[1]])
        for s_pos in np.linspace(0, self.s_max, 200):
            pos = P0[:2] + s_pos * n[:2]
            offset = self.r_beam * n_perp
            ax.plot([pos[0] - offset[0], pos[0] + offset[0]],
                    [pos[1] - offset[1], pos[1] + offset[1]],
                    'r-', linewidth=1.5, alpha=1)

        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_title('Flux Surfaces and NBI Beam (Top View)')
        ax.set_aspect('equal')
        ax.grid(True)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

    def todict(self):
        return {
            'enabled'    : self.enabled,
            's_max'      : self.s_max,
            'r_beam'     : self.r_beam,
            'P0'         : self.P0,
            'n'          : self.n,
            'energy_fractions' : self.energy_fractions,
            'Ti_beam'    : self.Ti_beam,
            'm_i_beam'   : self.m_i_beam,
            'R0'         : self.R0,
            'gaussian_profile' : self.gaussian_profile,
            'j_B'        : {
                't'       : self.j_B_t,
                'x'       : self.j_B_x,
                'tinterp' : self.j_B_tinterp,
            },
            'P_NBI'      : {
                't'       : self.P_NBI_t,
                'x'       : self.P_NBI_x,
                'tinterp' : self.P_NBI_tinterp,
            },
        }
