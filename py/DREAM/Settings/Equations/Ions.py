
import matplotlib.pyplot as plt
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.IonSpecies import IonSpecies, IONS_PRESCRIBED, IONIZATION_MODE_FLUID, IONIZATION_MODE_KINETIC, IONIZATION_MODE_KINETIC_APPROX_JAC, ION_OPACITY_MODE_TRANSPARENT
from . UnknownQuantity import UnknownQuantity

# Model to use for ion heat
IONS_T_I_NEGLECT = 1
IONS_T_I_INCLUDE = 2

class Ions(UnknownQuantity):
    

    def __init__(self, settings, ionization=IONIZATION_MODE_FLUID):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.ions = list()
        self.r    = None
        self.t    = None

        self.ionization = ionization
        self.typeTi = IONS_T_I_NEGLECT

    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, Z0=None, isotope=0, SPIMolarFraction=-1, opacity_mode=ION_OPACITY_MODE_TRANSPARENT, T=None, n=None, r=None, t=None, tritium=False):

        """
        Adds a new ion species to the plasma.

        :param str name:        Name by which the ion species will be referred to.
        :param int Z:           Ion charge number.
        :param int isotope:            Ion mass number.
        :param int iontype:     Method to use for evolving ions in time.
        :param int Z0:          Charge state to populate (used for populating exactly one charge state for the ion).
        :param n:               Ion density (can be either a scalar, 1D array or 2D array, depending on the other input parameters)
        :param float SPIMolarFraction: Molar fraction of the SPI injection (if any). A negative value means that this species is not part of the SPI injection 
        :param numpy.ndarray r: Radial grid on which the input density is defined.
        :param T:               Ion initial temperature (can be scalar for uniform temperature, otherwise 1D array matching `r` in size)
        :param numpy.ndarray r: Radial grid on which the input density and temperature is defined.
        :param numpy.ndarray t: Time grid on which the input density is defined.
        :param bool tritium:    If ``True``, the ion species is treated as Tritium.
        """
        if (self.r is not None) and (r is not None) and (np.any(self.r != r)):
            raise EquationException("The radial grid must be the same for all ion species.")
        if (self.t is not None) and (t is not None) and (np.any(self.t != t)):
            raise EquationException("The time grid must be the same for all ion species.")

        if T is not None:
            self.typeTi = IONS_T_I_INCLUDE

        ion = IonSpecies(settings=self.settings, name=name, Z=Z, ttype=iontype, Z0=Z0, isotope=isotope, SPIMolarFraction=SPIMolarFraction, opacity_mode=opacity_mode, T=T, n=n, r=r, t=t, interpr=self.r, interpt=None, tritium=tritium)

        self.ions.append(ion)

        self.r = ion.getR()
        if ion.getTime() is not None:
            self.t = ion.getTime()


    def getCharges(self):
        """
        Returns a list of the charges of the various ion species
        contained by this object.
        """
        return [ion.getZ() for ion in self.ions]

    def getIsotopes(self):
        """
        Returns a list of the isotopes of the various ion species
        contained by this object.
        """
        return [ion.getIsotope() for ion in self.ions]

    def getSPIMolarFraction(self):
        """
        Returns a list of the SPI molar fractions of the various ion species
        contained by this object.
        """
        return [ion.getSPIMolarFraction() for ion in self.ions]


    def getIon(self, i=None):
        """
        Returns the ion species with the specified index or name.

        :param i: Index or name of ion species to retrieve.
        """
        if type(i) == int: return self.ions[i]
        elif type(i) == str:
            for j in range(0, len(self.ions)):
                if self.ions[j].getName() == i:
                    return self.ions[j]

            raise EquationException("No ion with name '{}' has been defined.".format(i))
        else:
            raise EquationException("Invalid call to 'getIon()'.")


    def setIonization(self, ionization=IONIZATION_MODE_FLUID):
        """
        Sets which model to use for ionization.

        :param int ionization: Flag indicating which model to use for ionization.
        """
        self.ionization=ionization


    def getTritiumSpecies(self):
        """
        Returns a list of names of the ion species which are treated
        as Tritium.
        """
        trit = []
        for ion in self.ions:
            if ion.tritium:
                trit.append(ion.getName())

        return trit


    def getTypes(self):
        """
        Returns a list of ion types for the various ion species
        contained by this object.
        """
        return [ion.getType() for ion in self.ions]
        
    def getOpacityModes(self):
        """
        Returns a list of ion opacity modes for the various ion species
        contained by this object.
        """
        return [ion.getOpacityMode() for ion in self.ions]


    def setIonType(self, index, ttype):
        """
        Modifies the type of equation used for the specified ion species.

        :param index: Index or name of ion species to set type for.
        :param int ttype: Type of equation to use for evolving the ion species.
        """
        ion = self.getIon(index)

        # Note that the DREAM kernel only uses positive type indices.
        # The negative type indices are interface extensions which can
        # only be used with the 'initialize()' methods.
        if ttype <= 0:
            raise DREAMException("Trying to set invalid ion type for ion species '{}': {}.".format(ion.name, ttype))

        ion.ttype = ttype
    

    def fromdict(self, data):
        """
        Load settings from the specified dictionary.
        
        :param dict data: Dictionary containing all settings to load.
        """
        names        = data['names'].split(';')[:-1]
        Z            = data['Z']
        isotopes     = data['isotopes']
        types        = data['types']
        opacity_modes = data['opacity_modes']

        SPIMolarFraction = data['SPIMolarFraction']
        nZSPI = len(Z)-np.sum(SPIMolarFraction<0)
        if nZSPI>0:
            nShard = int(np.sum(SPIMolarFraction>=0)/nZSPI)
        else:
            nShard = 0

        if 'tritiumnames' in data:
            tritiumnames = data['tritiumnames'].split(';')[:-1]
        else:
            tritiumnames = []

        initial    = None
        prescribed = None
        initialTi  = None
        self.typeTi = IONS_T_I_NEGLECT
        if 'typeTi' in data:
            self.typeTi = int(data['typeTi'])
        if 'initial' in data:
            initial = data['initial']
        if 'prescribed' in data:
            prescribed = data['prescribed']
        if 'initialTi' in data:
            initialTi = data['initialTi']
        iidx, pidx, spiidx = 0, 0, 0
        for i in range(len(Z)):
            if types[i] == IONS_PRESCRIBED:
                n = prescribed['x'][pidx:(pidx+Z[i]+1)]
                r = prescribed['r']
                t = prescribed['t']
                pidx += Z[i]+1
            else:
                n = initial['x'][iidx:(iidx+Z[i]+1)]
                r = initial['r']
                t = None #initial['t']
                iidx += Z[i]+1
            if self.typeTi==IONS_T_I_INCLUDE and initialTi is not None:
                T = initialTi['x'][i]
            else: 
                T = None
            if SPIMolarFraction[spiidx]>=0:
                SPIMolarFractionSingleSpecies = SPIMolarFraction[spiidx:spiidx+nShard]
                spiidx+=nShard
            else:
                SPIMolarFractionSingleSpecies = SPIMolarFraction[spiidx]
                spiidx+=1
            tritium = (names[i] in tritiumnames)

            self.addIon(name=names[i], Z=Z[i], isotope=isotopes[i], SPIMolarFraction=SPIMolarFractionSingleSpecies, iontype=types[i], opacity_mode=opacity_modes[i], T=T, n=n, r=r, t=t, tritium=tritium)

        if 'ionization' in data:
            self.ionization = int(data['ionization'])

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this Ions object.
        """

        Z       = self.getCharges()
        itypes  = self.getTypes()
        iopacity_modes =self.getOpacityModes()
        isotopes     = self.getIsotopes()
        initial = None
        initialTi = None
        prescribed = None
        names   = ""

        tritiumnames = ""

        SPIMolarFraction = None

        for ion in self.ions:
            names += '{};'.format(ion.getName())

            if ion.tritium:
                tritiumnames += '{};'.format(ion.getName())

            if ion.getTime() is None:
                if initial is None:
                    initial = np.copy(ion.getDensity())
                else:
                    initial = np.concatenate((initial, ion.getDensity()))
            else:
                if prescribed is None:
                    prescribed = np.copy(ion.getDensity())
                else:
                    prescribed = np.concatenate((prescribed, ion.getDensity()))
            if initialTi is None:
                initialTi = np.copy(ion.getTemperature())
            else:
                initialTi = np.concatenate((initialTi, ion.getTemperature()))
                
            if SPIMolarFraction is None:
                SPIMolarFraction = np.copy(ion.getSPIMolarFraction())
            else:
                SPIMolarFraction = np.concatenate((SPIMolarFraction, ion.getSPIMolarFraction()))
                
        data = {
            'names': names,
            'Z': Z,
            'isotopes':isotopes,
            'SPIMolarFraction':SPIMolarFraction,
            'types': itypes,
            'opacity_modes':iopacity_modes
        }

        if len(tritiumnames) > 0:
            data['tritiumnames'] = tritiumnames

        if initial is not None:
            data['initial'] = {
                'r': self.r,
                'x': initial
            }

        if prescribed is not None:
            data['prescribed'] = {
                'r': self.r,
                't': self.t,
                'x': prescribed
            }

        data['initialTi'] = {
            'r': self.r,
            'x': initialTi
        }
        data['ionization'] = self.ionization
        data['typeTi'] = self.typeTi

        return data
            

    def verifySettings(self):
        """
        Verify that all settings are consistent.
        """
        # Make sure there are no double names
        for i in range(0, len(self.ions)):
            for j in range(0, len(self.ions)):
                if i == j: continue

                if self.ions[i].getName() == self.ions[j].getName():
                    raise EquationException("ions: More than one ion species is named '{}'.".format(self.ions[i].getName()))
            
            self.ions[i].verifySettings()
        
        if (self.ionization != IONIZATION_MODE_FLUID) and (self.ionization != IONIZATION_MODE_KINETIC) and (self.ionization != IONIZATION_MODE_KINETIC_APPROX_JAC):
            raise EquationException("ions: Invalid ionization mode: {}.".format(self.ionization))
 

    def getFreeElectronDensity(self, t=0):
        """
        Returns the plasma free electron density at the given time index, based
        on the prescribed/initialized ion densities.

        :param int t: Index of time for which to retrieve the free electron density.
        """
        n_free = np.zeros( self.r.shape )

        for ion in self.ions:
            for Z0 in range(1,ion.Z + 1):
                if len( ion.n.shape ) == 3:
                    n_free = n_free + Z0 * ion.n[Z0,t,:]
                elif len( ion.n.shape ) == 2:
                    n_free = n_free + Z0 * ion.n[Z0,:]
                
        return n_free, self.r


