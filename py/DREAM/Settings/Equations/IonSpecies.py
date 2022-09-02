#
# This class represents a single ion species (with multiple
# charge states). Each ion species provides a number of settings
# for DREAM:
#   
#   name     -- Label allowing user to identify ion species
#   Z        -- Ion charge number
#   type     -- Method to use for evolving the ion density
#               Is either 'prescribed' (time evolution explicitly
#               given in settings), 'equilibrium' (time evolution
#               governed by enforcing coronal equilibrium), or
#               'dynamic' (solve the ionization rate equation)
#   data
#     t      -- Time grid on which the density is given
#     r      -- Radial grid on which the density is given
#     n      -- Depending on 'type', this is either the prescribed
#               density (size (Z+1) x nt x nr) or the initial
#               density (size (Z+1) x nr), where 'nt' is the number
#               of time points that the density is given, and 'nr'
#               is the number of radial locations.
#

import matplotlib.pyplot as plt
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException


# Types in DREAM
IONS_PRESCRIBED = 1
IONS_EQUILIBRIUM = 2
IONS_DYNAMIC = 3

# Types which are extensions implemented in this interface
# (which are special cases of the DREAM types above)
IONS_DYNAMIC_NEUTRAL = -1
IONS_DYNAMIC_FULLY_IONIZED = -2
IONS_PRESCRIBED_NEUTRAL = -3
IONS_PRESCRIBED_FULLY_IONIZED = -4
IONS_EQUILIBRIUM_NEUTRAL = -5
IONS_EQUILIBRIUM_FULLY_IONIZED = -6

# Model to use for ionization
IONIZATION_MODE_FLUID = 1
IONIZATION_MODE_KINETIC = 2
IONIZATION_MODE_KINETIC_APPROX_JAC=3

ION_OPACITY_MODE_TRANSPARENT = 1
ION_OPACITY_MODE_GROUND_STATE_OPAQUE = 2

# Model to use for ion transport
ION_CHARGED_DIFFUSION_MODE_NONE = 1
ION_CHARGED_DIFFUSION_MODE_PRESCRIBED = 2

ION_NEUTRAL_DIFFUSION_MODE_NONE = 1
ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED = 2

ION_CHARGED_ADVECTION_MODE_NONE = 1
ION_CHARGED_ADVECTION_MODE_PRESCRIBED = 2

ION_NEUTRAL_ADVECTION_MODE_NONE = 1
ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED = 2

class IonSpecies:
    
    def __init__(self, settings, name, Z, ttype=0, Z0=None, isotope=0, SPIMolarFraction=-1.0, opacity_mode = ION_OPACITY_MODE_TRANSPARENT, 
        charged_diffusion_mode=ION_CHARGED_DIFFUSION_MODE_NONE, charged_prescribed_diffusion=None, rChargedPrescribedDiffusion=None, tChargedPrescribedDiffusion=None,
        neutral_diffusion_mode=ION_NEUTRAL_DIFFUSION_MODE_NONE, neutral_prescribed_diffusion=None, rNeutralPrescribedDiffusion=None, tNeutralPrescribedDiffusion=None,
        charged_advection_mode=ION_CHARGED_ADVECTION_MODE_NONE, charged_prescribed_advection=None, rChargedPrescribedAdvection=None, tChargedPrescribedAdvection=None,
        neutral_advection_mode=ION_NEUTRAL_ADVECTION_MODE_NONE, neutral_prescribed_advection=None, rNeutralPrescribedAdvection=None, tNeutralPrescribedAdvection=None,
        t_transp_expdecay_all_cs = None, t_transp_start_expdecay_all_cs = 0, diffusion_initial_all_cs = None, diffusion_final_all_cs = 0, advection_initial_all_cs = None, advection_final_all_cs = 0, r_expdecay_all_cs = None, t_expdecay_all_cs = None,        
        T=None, n=None, r=None, t=None, interpr=None, interpt=None, tritium=False, hydrogen=False):
        """
        Constructor.

        :param DREAMSettings settings: Parent DREAMSettings object.
        :param str name:               Name by which the ion species will be referred to.
        :param int Z:                  Ion charge number.
        :param int isotope:            Ion mass number.
        :param int ttype:              Method to use for evolving ions in time.
        :param int Z0:                 Charge state to populate with given density.
        :param float n:                Ion density (can be either a scalar, 1D array or 2D array, depending on the other input parameters)
        :param float SPIMolarFraction: Molar fraction of the SPI injection (if any). A negative value means that this species is not part of the SPI injection 
        :param T:                      Ion initial temperature (can be scalar for uniform temperature, otherwise 1D array matching `r` in size)
        :param numpy.ndarray r:        Radial grid on which the input density is defined.
        :param numpy.ndarray t:        Time grid on which the input density is defined.
        :param numpy.ndarray interpr:  Radial grid onto which ion densities should be interpolated.
        :param numpy.ndarray interpt:  Time grid onto which ion densities should be interpolated.
        :param bool tritium:           If ``True``, this ion species is treated as Tritium.
        :param bool hydrogen:          If ``True``, this ion species is treated as Hydrogen.
        """
        if ';' in name:
            raise EquationException("ion_species: '{}': Invalid character found in ion name: '{}'.".format(name, ';'))

        self.settings = settings
        self.name     = name
        self.Z        = int(Z)
        self.isotope  = int(isotope)
        self.ttype    = None
        self.tritium  = tritium
        self.hydrogen = hydrogen
        self.opacity_mode = opacity_mode
        self.charged_diffusion_mode = None
        self.neutral_diffusion_mode = None
        self.charged_advection_mode = None
        self.neutral_advection_mode = None

        self.setSPIMolarFraction(SPIMolarFraction)

        if self.tritium and self.hydrogen:
            raise EquationException("ion_species: '{}': Ion species indicated as both Tritium and Hydrogen simultaneously.")

        # Emit warning if 'T' is used as name but 'tritium = False',
        # as this may indicate a user error
        if name == 'T' and tritium == False:
            print("WARNING: Ion species with name 'T' added, but 'tritium = False'.")
        if name == 'H' and hydrogen == False:
            print("WARNING: Ion species with name 'H' added, but 'hydrogen = False'.")

        self.n = None
        self.r = None
        self.t = None
        if ttype == IONS_PRESCRIBED:
            if Z0 is not None:
                self.initialize_prescribed_charge_state(Z0=Z0, n=n, r=r, t=t, interpr=interpr, interpt=interpt)
            else:
                self.initialize_prescribed(n=n, r=r, t=t)
        elif ttype == IONS_DYNAMIC:
            if Z0 is not None:
                self.initialize_dynamic_charge_state(Z0=Z0, n=n, r=r, interpr=interpr)
            else:
                self.initialize_dynamic(n=n, r=r)
        elif ttype == IONS_EQUILIBRIUM:
            self.initialize_equilibrium(n=n, r=r, Z0=Z0)
        elif Z0 is not None:
            print("WARNING: Charge state Z0 given, but ion type is not simply 'prescribed', 'dynamic' or 'equilibrium'. Hence, Z0 is ignored.")
        
        # TYPES AVAILABLE ONLY IN THIS INTERFACE
        elif ttype == IONS_DYNAMIC_NEUTRAL:
            self.initialize_dynamic_neutral(n=n, r=r, interpr=interpr)
        elif ttype == IONS_DYNAMIC_FULLY_IONIZED:
            self.initialize_dynamic_fully_ionized(n=n, r=r, interpr=interpr)
        elif ttype == IONS_PRESCRIBED_NEUTRAL:
            self.initialize_prescribed_neutral(n=n, r=r, t=t, interpr=interpr, interpt=interpt)
        elif ttype == IONS_PRESCRIBED_FULLY_IONIZED:
            self.initialize_prescribed_fully_ionized(n=n, r=r, t=t, interpr=interpr, interpt=interpt)
        else:
            raise EquationException("ion_species: '{}': Unrecognized ion type: {}.".format(self.name, ttype))

        self.T = self.setTemperature(T)
        
        # Initialize diffusion
        self.charged_prescribed_diffusion = None
        self.rChargedPrescribedDiffusion = None
        self.tChargedPrescribedDiffusion = None
        if charged_diffusion_mode == ION_CHARGED_DIFFUSION_MODE_PRESCRIBED:
            # If an exponential decay of the transport coefficients are prescribed, 
            # set the precribed diffusion coefficients according to this, if nothing else is prescribed
            if charged_prescribed_diffusion is None and t_transp_expdecay_all_cs is not None:
                charged_prescribed_diffusion, rChargedPrescribedDiffusion, tChargedPrescribedDiffusion = self.calcTransportCoefficientExpdecayAllChargedStates(t_start = t_transp_start_expdecay_all_cs, t_exp = t_transp_expdecay_all_cs, c0 = diffusion_initial_all_cs, cf = diffusion_final_all_cs, r = r_expdecay_all_cs, t = t_expdecay_all_cs)
                
            self.initialize_charged_prescribed_diffusion(charged_prescribed_diffusion = charged_prescribed_diffusion, rChargedPrescribedDiffusion = rChargedPrescribedDiffusion,
                tChargedPrescribedDiffusion = tChargedPrescribedDiffusion, interpr=interpr, interpt=interpt)
        else:
            self.charged_diffusion_mode = charged_diffusion_mode
            
        self.neutral_prescribed_diffusion = None
        self.rNeutralPrescribedDiffusion = None
        self.tNeutralPrescribedDiffusion = None
        if neutral_diffusion_mode == ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED:
        
            # If an exponential decay of the transport coefficients are prescribed, 
            # set the precribed diffusion coefficients according to this, if nothing else is prescribed
            if neutral_prescribed_diffusion is None and t_transp_expdecay_all_cs is not None:
                neutral_prescribed_diffusion, rNeutralPrescribedDiffusion, tNeutralPrescribedDiffusion = self.calcTransportCoefficientExpdecaySingleChargeState(t_start = t_transp_start_expdecay_all_cs, t_exp = t_transp_expdecay_all_cs, c0 = diffusion_initial_all_cs, cf = diffusion_final_all_cs, r = r_expdecay_all_cs, t = t_expdecay_all_cs)
                
            self.initialize_neutral_prescribed_diffusion(neutral_prescribed_diffusion = neutral_prescribed_diffusion, rNeutralPrescribedDiffusion = rNeutralPrescribedDiffusion,
                tNeutralPrescribedDiffusion = tNeutralPrescribedDiffusion, interpr=interpr, interpt=interpt)
        else:
            self.neutral_diffusion_mode = neutral_diffusion_mode

        # Initialize advection
        self.charged_prescribed_advection = None
        self.rChargedPrescribedAdvection = None
        self.tChargedPrescribedAdvection = None
        if charged_advection_mode == ION_CHARGED_ADVECTION_MODE_PRESCRIBED:
        
            # If an exponential decay of the transport coefficients are prescribed, 
            # set the precribed advection coefficients according to this, if nothing else is prescribed
            if charged_prescribed_advection is None and t_transp_expdecay_all_cs is not None:
                charged_prescribed_advection, rChargedPrescribedAdvection, tChargedPrescribedAdvection = self.calcTransportCoefficientExpdecayAllChargedStates(t_start = t_transp_start_expdecay_all_cs, t_exp = t_transp_expdecay_all_cs, c0 = advection_initial_all_cs, cf = advection_final_all_cs, r = r_expdecay_all_cs, t = t_expdecay_all_cs)
                
            self.initialize_charged_prescribed_advection(charged_prescribed_advection = charged_prescribed_advection, rChargedPrescribedAdvection = rChargedPrescribedAdvection,
                tChargedPrescribedAdvection = tChargedPrescribedAdvection, interpr=interpr, interpt=interpt)
        else:
            self.charged_advection_mode = charged_advection_mode
            
        self.neutral_prescribed_advection = None
        self.rNeutralPrescribedAdvection = None
        self.tNeutralPrescribedAdvection = None
        if neutral_advection_mode == ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED:
        
            # If an exponential decay of the transport coefficients are prescribed, 
            # set the precribed advection coefficients according to this, if nothing else is prescribed
            if neutral_prescribed_advection is None and t_transp_expdecay_all_cs is not None:
                neutral_prescribed_advection, rNeutralPrescribedAdvection, tNeutralPrescribedAdvection = self.calcTransportCoefficientExpdecaySingleChargeState(t_start = t_transp_start_expdecay_all_cs, t_exp = t_transp_expdecay_all_cs, c0 = advection_initial_all_cs, cf = advection_final_all_cs, r = r_expdecay_all_cs, t = t_expdecay_all_cs)
                
            self.initialize_neutral_prescribed_advection(neutral_prescribed_advection = neutral_prescribed_advection, rNeutralPrescribedAdvection = rNeutralPrescribedAdvection,
                tNeutralPrescribedAdvection = tNeutralPrescribedAdvection, interpr=interpr, interpt=interpt)
        else:
            self.neutral_advection_mode = neutral_advection_mode

    def setTemperature(self, T):
        """
        Sets the ion temperature from an input value `T`. 
        For scalar T, sets a uniform radial profile,
        otherwise requires the T profile to be given on the 
        `r` grid which is provided to the IonSpecies constructor.
        """
        if type(T) == list:
            T = np.array(T)
        if T is None:
            T = np.zeros((1,np.size(self.r)))
        elif np.isscalar(T):
            T = np.ones((1, np.size(self.r)))*T
        elif np.ndim(T)==1:  
            T = T[None,:]
        elif T.shape[1] != np.size(self.r):
             raise EquationException("ion_species: '{}': Invalid dimensions of initial ion temperature T: {}x{}. Expected {}x{}."
                .format(self.name, T.shape[0], T.shape[1], 1, np.size(self.r)))        
        return T
        
    def shiftTimeTranspCoeffs(self, tShift):
        if self.tChargedPrescribedDiffusion is not None:
            self.tChargedPrescribedDiffusion = self.tChargedPrescribedDiffusion - tShift
        if self.tNeutralPrescribedDiffusion is not None:    
            self.tNeutralPrescribedDiffusion = self.tNeutralPrescribedDiffusion - tShift
        if self.tChargedPrescribedAdvection is not None:
            self.tChargedPrescribedAdvection = self.tChargedPrescribedAdvection - tShift
        if self.tNeutralPrescribedAdvection is not None:
            self.tNeutralPrescribedAdvection = self.tNeutralPrescribedAdvection - tShift

    def getDensity(self):
        """
        Returns the prescribed density array for this ion species.
        """
        return self.n
        
        
        
    # Getters for diffusion-related quantities    
    def getChargedPrescribedDiffusion(self):
        """
        Returns the prescribed charged diffusion coefficient array for this ion species.
        """
        return self.charged_prescribed_diffusion
        
    def getRChargedPrescribedDiffusion(self):
        """
        Returns the radial grid for the prescribed charged diffusion coefficient array for this ion species.
        """
        return self.rChargedPrescribedDiffusion
        
    def getTChargedPrescribedDiffusion(self):
        """
        Returns the time grid for the prescribed charged diffusion coefficient array for this ion species.
        """
        return self.tChargedPrescribedDiffusion
        
    def getRNeutralPrescribedDiffusion(self):
        """
        Returns the radial grid for the prescribed neutral diffusion coefficient array for this ion species.
        """
        return self.rNeutralPrescribedDiffusion
        
    def getTNeutralPrescribedDiffusion(self):
        """
        Returns the time grid for the prescribed neutral diffusion coefficient array for this ion species.
        """
        return self.tNeutralPrescribedDiffusion
        
    def getNeutralPrescribedDiffusion(self):
        """
        Returns the prescribed neutral diffusion coefficient array for this ion species.
        """
        return self.neutral_prescribed_diffusion
        
        
        
    # Getters for advection-related quantities    
    def getChargedPrescribedAdvection(self):
        """
        Returns the prescribed charged advection coefficient array for this ion species.
        """
        return self.charged_prescribed_advection
        
    def getRChargedPrescribedAdvection(self):
        """
        Returns the radial grid for the prescribed charged advection coefficient array for this ion species.
        """
        return self.rChargedPrescribedAdvection
        
    def getTChargedPrescribedAdvection(self):
        """
        Returns the time grid for the prescribed charged advection coefficient array for this ion species.
        """
        return self.tChargedPrescribedAdvection
        
    def getRNeutralPrescribedAdvection(self):
        """
        Returns the radial grid for the prescribed neutral advection coefficient array for this ion species.
        """
        return self.rNeutralPrescribedAdvection
        
    def getTNeutralPrescribedAdvection(self):
        """
        Returns the time grid for the prescribed neutral advection coefficient array for this ion species.
        """
        return self.tNeutralPrescribedAdvection
        
    def getNeutralPrescribedAdvection(self):
        """
        Returns the prescribed neutral advection coefficient array for this ion species.
        """
        return self.neutral_prescribed_advection
        
        

    def getName(self):
        """
        Returns the name of this ion species.
        """
        return self.name


    def getR(self):
        """
        Returns the radial grid on which the ion densities are defined.
        """
        return self.r


    def getTime(self):
        """
        Returns the time grid on which the ion densities are defined.
        """
        return self.t


    def getType(self):
        """
        Returns the type of equation to use for evolving the ion densities
        for this species.
        """
        return self.ttype
        
    def getOpacityMode(self):
        """
        Returns the opacity mode to use for evolving the ion densities
        for this species.
        """
        return self.opacity_mode
        
    def getChargedDiffusionMode(self):
        """
        Returns the charged diffusion mode to use for evolving the ion densities
        for this species.
        """
        return self.charged_diffusion_mode
        
    def getNeutralDiffusionMode(self):
        """
        Returns the neutral diffusion mode to use for evolving the ion densities
        for this species.
        """
        return self.neutral_diffusion_mode
        
    def getChargedAdvectionMode(self):
        """
        Returns the charged advection mode to use for evolving the ion densities
        for this species.
        """
        return self.charged_advection_mode
        
    def getNeutralAdvectionMode(self):
        """
        Returns the neutral advection mode to use for evolving the ion densities
        for this species.
        """
        return self.neutral_advection_mode

    def getTemperature(self):
        """
        Returns the initial temperature array to use for evolving
        the ion heat of this species 
        """
        return self.T

    def getZ(self):
        """
        Returns the atomic charge for this ion species.
        """
        return self.Z


    def getIsotope(self): return self.isotope


    def getSPIMolarFraction(self): return self.SPIMolarFraction
    
    def setSPIMolarFraction(self, SPIMolarFraction):
        if np.isscalar(SPIMolarFraction):
            self.SPIMolarFraction = np.array([SPIMolarFraction])
        else:
            self.SPIMolarFraction = SPIMolarFraction


    def isHydrogen(self):
        """
        Returns ``True`` if this ion species is a hydrogen species.
        """
        return self.hydrogen


    def isTritium(self):
        """
        Returns ``True`` if this ion species is a tritium species.
        """
        return self.tritium


    def initialize_prescribed(self, n=None, r=None, t=None):
        """
        Prescribes the evolution for this ion species.
        """
        self.ttype = IONS_PRESCRIBED
        if n is None:
            raise EquationException("ion_species: '{}': Input density must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(n) == list:
            n = np.array(n)

        # Scalar (assume density constant in spacetime)
        #if type(n) == float or (type(n) == np.ndarray and n.size == 1):
        if np.isscalar(n):
            self.t = np.array([0])
            self.r = np.array([0,1])
            self.n = np.ones((self.Z+1,1,2)) * n
            return
        if r is None:
            raise EquationException("ion_species: '{}': Non-scalar density prescribed, but no radial coordinates given.".format(self.name))

        # Radial profile (assume fully ionized)
        if len(n.shape) == 1:
            raise EquationException("ion_species: '{}': Prescribed density data has only one dimension.".format(self.name))
        # Radial profiles of charge states
        elif len(n.shape) == 2:
            raise EquationException("ion_species: '{}': Prescribed density data has only two dimensions.".format(self.name))
        # Full time evolution of radial profiles of charge states
        elif len(n.shape) == 3:
            if t is None:
                raise EquationException("ion_species: '{}': 3D ion density prescribed, but no time coordinates given.".format(self.name))

            if self.Z+1 != n.shape[0] or t.size != n.shape[1] or r.size != n.shape[2]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed density: {}x{}x{}. Expected {}x{}x{}"
                    .format(self.name, n.shape[0], n.shape[1], n.shape[2], self.Z+1, t.size, r.size))
            self.t = t
            self.r = r
            self.n = n
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of prescribed density: {}.".format(self.name, n.shape))


    def initialize_dynamic(self, n=None, r=None):
        """
        Evolve ions according to the ion rate equation in DREAM.
        """
        self.ttype = IONS_DYNAMIC

        if n is None:
            raise EquationException("ion_species: '{}': Input density must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(n) == list:
            n = np.array(n)

        # Scalar (assume density constant in spacetime)
        if type(n) == float or (type(n) == np.ndarray and n.size == 1):
            raise EquationException("ion_species: '{}': Initial density must be two dimensional (charge states x radius).".format(self.name))

        if r is None:
            raise EquationException("ion_species: '{}': Non-scalar initial ion density prescribed, but no radial coordinates given.".format(self.name))

        # Radial profiles for all charge states 
        if len(n.shape) == 2:
            if self.Z+1 != n.shape[0] or r.size != n.shape[1]:
                raise EquationException("ion_species: '{}': Invalid dimensions of initial ion density: {}x{}. Expected {}x{}."
                    .format(self.name, n.shape[0], n.shape[1], self.Z+1, r.size))

            self.t = None
            self.r = r
            self.n = n
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of initial density: {}.".format(n.shape).format(self.name))


    def initialize_equilibrium(self, n=None, r=None, interpr=None):
        """
        Evolve ions according to the equilibrium equation in DREAM.
        """
        self.ttype = IONS_EQUILIBRIUM

        if n is None:
            raise EquationException("ion_species: '{}': Input density must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(n) == list:
            n = np.array(n)

        # Scalar (assume density constant in radius)
        if type(n) == float or (type(n) == np.ndarray and n.size == 1):
            r = interpr if interpr is not None else np.array([0])
            N = np.zeros((self.Z+1,r.size))

            # For the equilibrium, it doesn't matter which charge state we
            # put the particles in. They will be placed in the correct state
            # after the first time step (=> make all particles fully ionized
            # so that nfree > 0)
            N[self.Z,:] = n
            n = N
        elif r is None:
            raise EquationException("ion_species: '{}': Non-scalar initial ion density prescribed, but no radial coordinates given.".format(self.name))

        # Radial profiles for all charge states 
        if len(n.shape) == 2:
            if self.Z+1 != n.shape[0] or r.size != n.shape[1]:
                raise EquationException("ion_species: '{}': Invalid dimensions of initial ion density: {}x{}. Expected {}x{}."
                    .format(self.name, n.shape[0], n.shape[1], self.Z+1, r.size))

            self.t = None
            self.r = r
            self.n = n
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of initial density: {}.".format(self.name, n.shape))


    def initialize_dynamic_neutral(self, n=None, r=None, interpr=None):
        """
        Evolve the ions dynamically, initializing them all as neutrals.
        """
        self.initialize_dynamic_charge_state(0, n=n, r=r, interpr=interpr)


    def initialize_dynamic_fully_ionized(self, n=None, r=None, interpr=None):
        """
        Evolve the ions dynamically, initializing them all as fully ionized.
        """
        self.initialize_dynamic_charge_state(self.Z, n=n, r=r, interpr=interpr)


    def initialize_dynamic_charge_state(self, Z0, n=None, r=None, interpr=None):
        """
        Evolve the ions dynamically, initializing them all to reside in the specified charge state Z0.
        """
        if Z0 > self.Z or Z0 < 0:
            raise EquationException("ion_species: '{}': Invalid charge state specified: {}. Ion has charge Z = {}.".format(self.name, Z0, self.Z))

        if n is None:
            raise EquationException("ion_species: '{}': Input density must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(n) == list:
            n = np.array(n)

        # Scalar (assume density constant in spacetime)
        if type(n) == float or np.isscalar(n) or (type(n) == np.ndarray and n.size == 1):
            r = interpr if interpr is not None else np.array([0])
            N = np.zeros((self.Z+1,r.size))
            N[Z0,:] = n

            self.initialize_dynamic(n=N, r=r)
            return

        if r is None:
            raise EquationException("ion_species: '{}': Non-scalar density prescribed, but no radial coordinates given.".format(self.name))

        # Radial profile
        if len(n.shape) == 1:
            if r.size != n.size:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed density: {}. Expected {}."
                    .format(self.name, n.shape[0], r.size))
                
            N = np.zeros((self.Z+1, r.size))
            N[Z0,:] = n
            self.initialize_dynamic(n=N, r=r)
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of prescribed density: {}.".format(self.name, n.shape))


    def initialize_prescribed_neutral(self, n=None, r=None, t=None, interpr=None, interpt=None):
        """
        Prescribe the ions to be neutral.
        """
        self.initialize_prescribed_charge_state(0, n=n, r=r, t=t, interpr=interpr, interpt=interpt)


    def initialize_prescribed_fully_ionized(self, n=None, r=None, t=None, interpr=None, interpt=None):
        """
        Prescribe the ions to be fully ionized.
        """
        self.initialize_prescribed_charge_state(self.Z, n=n, r=r, t=t, interpr=interpr, interpt=interpt)


    def initialize_prescribed_charge_state(self, Z0, n=None, r=None, t=None, interpr=None, interpt=None):
        """
        Prescribe the ions to all be situated in the specified charge state Z0.
        """
        if Z0 > self.Z or Z0 < 0:
            raise EquationException("ion_species: '{}': Invalid charge state specified: {}. Ion has charge Z = {}.".format(self.name, Z0, self.Z))

        if n is None:
            raise EquationException("ion_species: '{}': Input density must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(n) == list:
            n = np.array(n)

        # Scalar (assume density constant in spacetime)
        #if type(n) == float or (type(n) == np.ndarray and n.size == 1):
        if np.isscalar(n):
            t = interpt if interpt is not None else np.array([0])
            r = interpr if interpr is not None else np.array([0])
            N = np.zeros((self.Z+1,t.size,r.size))
            N[Z0,0,:] = n

            self.initialize_prescribed(n=N, t=t, r=r)
            return

        if r is None:
            raise EquationException("ion_species: '{}': Non-scalar density prescribed, but no radial coordinates given.".format(self.name))

        # Radial profile
        if len(n.shape) == 1:
            if r.size != n.size:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed density: {}. Expected {}."
                    .format(self.name, n.shape[0], r.size))
                
            t = interpt if interpt is not None else np.array([0])
            n = np.reshape(n, (t.size,r.size))

        # Radial + temporal profile
        if len(n.shape) == 2:
            if t is None:
                raise EquationException("ion_species: '{}': 2D ion density prescribed, but no time coordinates given.".format(self.name))

            if t.size != n.shape[0] or r.size != n.shape[1]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed density: {}x{}. Expected {}x{}."
                    .format(self.name, n.shape[0], n.shape[1], t.size, r.size))

            N = np.zeros((self.Z+1, t.size, r.size))
            N[Z0,:,:] = n

            self.initialize_prescribed(n=N, t=t, r=r)
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of prescribed density: {}.".format(self.name, n.shape))
            
            
            
            
    def initialize_charged_prescribed_diffusion(self, charged_prescribed_diffusion=None, rChargedPrescribedDiffusion=None, tChargedPrescribedDiffusion=None, interpr=None, interpt=None):
        """
        Prescribes the evolution of the charged diffusion coefficients for this ion species.
        """
        self.charged_diffusion_mode = ION_CHARGED_DIFFUSION_MODE_PRESCRIBED
        if charged_prescribed_diffusion is None:
            raise EquationException("ion_species: '{}': Prescribed charged diffusion coefficients must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(charged_prescribed_diffusion) == list:
            charged_prescribed_diffusion = np.array(charged_prescribed_diffusion)

        # Scalar (assume density constant in spacetime)
        if np.isscalar(charged_prescribed_diffusion):
            self.tChargedPrescribedDiffusion = np.array([0])
            self.rChargedPrescribedDiffusion = np.array([0,1])
            self.charged_prescribed_diffusion = np.ones((self.Z,1,2)) * charged_prescribed_diffusion
            return
        if rChargedPrescribedDiffusion is None:
            raise EquationException("ion_species: '{}': Non-scalar density prescribed, but no radial coordinates given.".format(self.name))

        if len(charged_prescribed_diffusion.shape) == 1:
            raise EquationException("ion_species: '{}': Prescribed charged diffusion coefficient data has only one dimension.".format(self.name))
        elif len(charged_prescribed_diffusion.shape) == 2:
            raise EquationException("ion_species: '{}': Prescribed charged diffusion coefficient data has only two dimensions.".format(self.name))
        # Full time evolution of radial profiles of charge states
        elif len(charged_prescribed_diffusion.shape) == 3:
            if tChargedPrescribedDiffusion is None:
                raise EquationException("ion_species: '{}': 3D charged diffusion coefficient prescribed, but no time coordinates given.".format(self.name))

            if self.Z != charged_prescribed_diffusion.shape[0] or tChargedPrescribedDiffusion.size != charged_prescribed_diffusion.shape[1] or rChargedPrescribedDiffusion.size != charged_prescribed_diffusion.shape[2]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed charged diffusion coefficient: {}x{}x{}. Expected {}x{}x{}"
                    .format(self.name, charged_prescribed_diffusion.shape.shape[0], charged_prescribed_diffusion.shape.shape[1], charged_prescribed_diffusion.shape.shape[2], self.Z, tChargedPrescribedDiffusion.size, rChargedPrescribedDiffusion.size))
            self.tChargedPrescribedDiffusion = tChargedPrescribedDiffusion
            self.rChargedPrescribedDiffusion = rChargedPrescribedDiffusion
            self.charged_prescribed_diffusion = charged_prescribed_diffusion
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of prescribed charged diffusion coefficient: {}.".format(self.name, charged_prescribed_diffusion.shape))

    def initialize_neutral_prescribed_diffusion(self, neutral_prescribed_diffusion=None, rNeutralPrescribedDiffusion=None, tNeutralPrescribedDiffusion=None, interpr=None, interpt=None):
        """
        Prescribes the evolution of the neutral diffusion coefficients for this ion species.
        """
        self.neutral_diffusion_mode = ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED
        if neutral_prescribed_diffusion is None:
            raise EquationException("ion_species: '{}': Prescribed neutral diffusion coefficients must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(neutral_prescribed_diffusion) == list:
            neutral_prescribed_diffusion = np.array(neutral_prescribed_diffusion)

        # Scalar (assume density constant in spacetime)
        if np.isscalar(neutral_prescribed_diffusion):
            self.tNeutralPrescribedDiffusion = np.array([0])
            self.rNeutralPrescribedDiffusion = np.array([0,1])
            self.neutral_prescribed_diffusion = np.ones((1,1,2)) * neutral_prescribed_diffusion
            return
        if rNeutralPrescribedDiffusion is None:
            raise EquationException("ion_species: '{}': Non-scalar density prescribed, but no radial coordinates given.".format(self.name))

        if len(neutral_prescribed_diffusion.shape) == 1:
            raise EquationException("ion_species: '{}': Prescribed neutral diffusion coefficient data has only one dimension.".format(self.name))
        # As there is only one neutral charge state for a single species, all information needed here can actually be provided in a 2D array
        elif len(neutral_prescribed_diffusion.shape) == 2:
            if tNeutralPrescribedDiffusion is None:
                raise EquationException("ion_species: '{}': 2D neutral diffusion coefficient prescribed, but no time coordinates given.".format(self.name))

            if tNeutralPrescribedDiffusion.size != neutral_prescribed_diffusion.shape[0] or rNeutralPrescribedDiffusion.size != neutral_prescribed_diffusion.shape[1]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed neutral diffusion coefficient: {}x{}. Expected {}x{}"
                    .format(self.name, neutral_prescribed_diffusion.shape.shape[0], neutral_prescribed_diffusion.shape.shape[1], tNeutralPrescribedDiffusion.size, rNeutralPrescribedDiffusion.size))
            self.tNeutralPrescribedDiffusion = tNeutralPrescribedDiffusion
            self.rNeutralPrescribedDiffusion = rNeutralPrescribedDiffusion
            self.neutral_prescribed_diffusion = neutral_prescribed_diffusion.reshape((1,neutral_prescribed_diffusion.shape[0],neutral_prescribed_diffusion.shape[1]))
        # Full time evolution of radial profiles of charge states
        elif len(neutral_prescribed_diffusion.shape) == 3:
            if tNeutralPrescribedDiffusion is None:
                raise EquationException("ion_species: '{}': 3D neutral diffusion coefficient prescribed, but no time coordinates given.".format(self.name))

            if neutral_prescribed_diffusion.shape[0] != 1 or tNeutralPrescribedDiffusion.size != neutral_prescribed_diffusion.shape[1] or rNeutralPrescribedDiffusion.size != neutral_prescribed_diffusion.shape[2]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed neutral diffusion coefficient: {}x{}x{}. Expected {}x{}x{}"
                    .format(self.name, neutral_prescribed_diffusion.shape[0], neutral_prescribed_diffusion.shape[1], neutral_prescribed_diffusion.shape[2], self.Z, tNeutralPrescribedDiffusion.size, rNeutralPrescribedDiffusion.size))
            self.tNeutralPrescribedDiffusion = tNeutralPrescribedDiffusion
            self.rNeutralPrescribedDiffusion = rNeutralPrescribedDiffusion
            self.neutral_prescribed_diffusion = neutral_prescribed_diffusion
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of prescribed neutral diffusion coefficient: {}.".format(self.name, neutral_prescribed_diffusion.shape))
   
   
   
   
    def initialize_charged_prescribed_advection(self, charged_prescribed_advection=None, rChargedPrescribedAdvection=None, tChargedPrescribedAdvection=None, interpr=None, interpt=None):
        """
        Prescribes the evolution of the charged advection coefficients for this ion species.
        """
        self.charged_advection_mode = ION_CHARGED_ADVECTION_MODE_PRESCRIBED
        if charged_prescribed_advection is None:
            raise EquationException("ion_species: '{}': Prescribed charged advection coefficients must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(charged_prescribed_advection) == list:
            charged_prescribed_advection = np.array(charged_prescribed_advection)

        # Scalar (assume density constant in spacetime)
        if np.isscalar(charged_prescribed_advection):
            self.tChargedPrescribedAdvection = np.array([0])
            self.rChargedPrescribedAdvection = np.array([0,1])
            self.charged_prescribed_advection = np.ones((self.Z,1,2)) * charged_prescribed_advection
            return
        if rChargedPrescribedAdvection is None:
            raise EquationException("ion_species: '{}': Non-scalar density prescribed, but no radial coordinates given.".format(self.name))

        if len(charged_prescribed_advection.shape) == 1:
            raise EquationException("ion_species: '{}': Prescribed charged advection coefficient data has only one dimension.".format(self.name))
        elif len(charged_prescribed_advection.shape) == 2:
            raise EquationException("ion_species: '{}': Prescribed charged advection coefficient data has only two dimensions.".format(self.name))
        # Full time evolution of radial profiles of charge states
        elif len(charged_prescribed_advection.shape) == 3:
            if tChargedPrescribedAdvection is None:
                raise EquationException("ion_species: '{}': 3D charged advection coefficient prescribed, but no time coordinates given.".format(self.name))

            if self.Z != charged_prescribed_advection.shape[0] or tChargedPrescribedAdvection.size != charged_prescribed_advection.shape[1] or rChargedPrescribedAdvection.size != charged_prescribed_advection.shape[2]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed charged advection coefficient: {}x{}x{}. Expected {}x{}x{}"
                    .format(self.name, charged_prescribed_advection.shape.shape[0], charged_prescribed_advection.shape.shape[1], charged_prescribed_advection.shape.shape[2], self.Z, tChargedPrescribedAdvection.size, rChargedPrescribedAdvection.size))
            self.tChargedPrescribedAdvection = tChargedPrescribedAdvection
            self.rChargedPrescribedAdvection = rChargedPrescribedAdvection
            self.charged_prescribed_advection = charged_prescribed_advection
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of prescribed charged advection coefficient: {}.".format(self.name, charged_prescribed_advection.shape))

    def initialize_neutral_prescribed_advection(self, neutral_prescribed_advection=None, rNeutralPrescribedAdvection=None, tNeutralPrescribedAdvection=None, interpr=None, interpt=None):
        """
        Prescribes the evolution of the neutral advection coefficients for this ion species.
        """
        self.neutral_advection_mode = ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED
        if neutral_prescribed_advection is None:
            raise EquationException("ion_species: '{}': Prescribed neutral advection coefficients must not be 'None'.".format(self.name))

        # Convert lists to NumPy arrays
        if type(neutral_prescribed_advection) == list:
            neutral_prescribed_advection = np.array(neutral_prescribed_advection)

        # Scalar (assume density constant in spacetime)
        if np.isscalar(neutral_prescribed_advection):
            self.tNeutralPrescribedAdvection = np.array([0])
            self.rNeutralPrescribedAdvection = np.array([0,1])
            self.neutral_prescribed_advection = np.ones((1,1,2)) * neutral_prescribed_advection
            return
        if rNeutralPrescribedAdvection is None:
            raise EquationException("ion_species: '{}': Non-scalar density prescribed, but no radial coordinates given.".format(self.name))

        if len(neutral_prescribed_advection.shape) == 1:
            raise EquationException("ion_species: '{}': Prescribed neutral advection coefficient data has only one dimension.".format(self.name))
        # As there is only one neutral charge state for a single species, all information needed here can actually be provided in a 2D array
        elif len(neutral_prescribed_advection.shape) == 2:
            if tNeutralPrescribedAdvection is None:
                raise EquationException("ion_species: '{}': 2D neutral advection coefficient prescribed, but no time coordinates given.".format(self.name))

            if tNeutralPrescribedAdvection.size != neutral_prescribed_advection.shape[0] or rNeutralPrescribedAdvection.size != neutral_prescribed_advection.shape[1]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed neutral advection coefficient: {}x{}. Expected {}x{}"
                    .format(self.name, neutral_prescribed_advection.shape.shape[0], neutral_prescribed_advection.shape.shape[1], tNeutralPrescribedAdvection.size, rNeutralPrescribedAdvection.size))
            self.tNeutralPrescribedAdvection = tNeutralPrescribedAdvection
            self.rNeutralPrescribedAdvection = rNeutralPrescribedAdvection
            self.neutral_prescribed_advection = neutral_prescribed_advection.reshape((1,neutral_prescribed_advection.shape[0],neutral_prescribed_advection.shape[1]))
        # Full time evolution of radial profiles of charge states
        elif len(neutral_prescribed_advection.shape) == 3:
            if tNeutralPrescribedAdvection is None:
                raise EquationException("ion_species: '{}': 3D neutral advection coefficient prescribed, but no time coordinates given.".format(self.name))

            if neutral_prescribed_advection.shape[0] != 1 or tNeutralPrescribedAdvection.size != neutral_prescribed_advection.shape[1] or rNeutralPrescribedAdvection.size != neutral_prescribed_advection.shape[2]:
                raise EquationException("ion_species: '{}': Invalid dimensions of prescribed neutral advection coefficient: {}x{}x{}. Expected {}x{}x{}"
                    .format(self.name, neutral_prescribed_advection.shape.shape[0], neutral_prescribed_advection.shape.shape[1], neutral_prescribed_advection.shape.shape[2], self.Z, tNeutralPrescribedAdvection.size, rNeutralPrescribedAdvection.size))
            self.tNeutralPrescribedAdvection = tNeutralPrescribedAdvection
            self.rNeutralPrescribedAdvection = rNeutralPrescribedAdvection
            self.neutral_prescribed_advection = neutral_prescribed_advection
        else:
            raise EquationException("ion_species: '{}': Unrecognized shape of prescribed neutral advection coefficient: {}.".format(self.name, neutral_prescribed_advection.shape))


    def calcTransportCoefficientExpdecaySingleChargeState(self, t_exp, c0, cf = 0, t_start = 0, r = None, t = None):
        if t is None:
            t = np.linspace(0,t_start+10*t_exp).reshape(-1,1)
        if r is None:
            r = np.linspace(0,self.settings.radialgrid.a)
        if np.isscalar(c0):
            Nr = len(r)
            c0 = c0*np.ones((1,Nr))
            
        if np.isscalar(cf):
            Nr = len(r)
            cf = cf*np.ones((1,Nr))     
                    
        c_single_charge_state = (cf + np.exp(-(t-t_start)/t_exp)*(c0-cf))*(t>t_start)
        
        return c_single_charge_state, r.flatten(), t.flatten()

    def calcTransportCoefficientExpdecayAllChargedStates(self, t_exp, c0, cf = 0, t_start = 0, r = None, t = None):
        c_single_charge_state, r, t = self.calcTransportCoefficientExpdecaySingleChargeState(t_exp, c0, cf, t_start, r, t)
        cCharged = np.zeros((self.Z,len(t),len(c_single_charge_state)))
        for i in range(self.Z):
            cCharged[i,:,:]=c_single_charge_state
        
        return cCharged, r, t

    
    def verifySettings(self):
        """
        Verify that the settings of this ion species are correctly set.
        """
        if self.Z < 1:
            raise EquationException("ion_species: '{}': Invalid atomic charge: {}.".format(self.name, self.Z))

        if self.hydrogen or self.tritium:
            if self.Z != 1:
                raise EquationException(f"ion_species: '{self.name}': Ion indicated as Hydrogen/Tritium, but charge Z = {self.Z}.")

        if self.ttype == IONS_PRESCRIBED:
            if self.t.ndim != 1:
                raise EquationException("ion_species: '{}': The time vector must be 1D.".format(self.name))
            elif self.r.ndim != 1:
                raise EquationException("ion_species: '{}': The time vector must be 1D.".format(self.name))
            elif self.n is None or (self.n.shape != (self.Z+1, self.t.size, self.r.size)):
                raise EquationException("ion_species: '{}': Invalid dimensions for input density: {}x{}x{}. Expected {}x{}x{}."
                    .format(self.name, self.n.shape[0], self.n.shape[1], self.n.shape[2], self.Z+1, self.t.size, self.r.size))
        elif self.ttype == IONS_EQUILIBRIUM or self.ttype == IONS_DYNAMIC:
            if (self.r is None) or (self.r.ndim != 1):
                raise EquationException("ion_species: '{}': The time vector must be 1D.".format(self.name))
            elif (self.n is None) or (self.n.shape != (self.Z+1, self.r.size)):
                raise EquationException("ion_species: '{}': Invalid dimensions for input density: {}x{}. Expected {}x{}."
                    .format(self.name, self.n.shape[0], self.n.shape[1], self.Z+1, self.r.size))


