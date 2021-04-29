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

class IonSpecies:
    
    def __init__(self, settings, name, Z, ttype=0, Z0=None, isotope=0, SPIMolarFraction=-1.0, opacity_mode = ION_OPACITY_MODE_TRANSPARENT, T=None, n=None, r=None, t=None, interpr=None, interpt=None, tritium=False):
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
        """
        if ';' in name:
            raise EquationException("ion_species: '{}': Invalid character found in ion name: '{}'.".format(name, ';'))

        self.settings = settings
        self.name     = name
        self.Z        = int(Z)
        self.isotope  = int(isotope)
        self.ttype    = None
        self.tritium  = tritium
        self.opacity_mode = opacity_mode

        if np.isscalar(SPIMolarFraction):
            self.SPIMolarFraction = np.array([SPIMolarFraction])
        else:
            self.SPIMolarFraction = SPIMolarFraction

        # Emit warning if 'T' is used as name but 'tritium = False',
        # as this may indicate a user error
        if name == 'T' and tritium == False:
            print("WARNING: Ion species with name 'T' added, but 'tritium = False'.")
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

    def getDensity(self):
        """
        Returns the prescribed density array for this ion species.
        """
        return self.n

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

    
    def verifySettings(self):
        """
        Verify that the settings of this ion species are correctly set.
        """
        if self.Z < 1:
            raise EquationException("ion_species: '{}': Invalid atomic charge: {}.".format(self.Z))

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


