
import matplotlib.pyplot as plt
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.IonSpecies import IonSpecies, IONS_PRESCRIBED
from . UnknownQuantity import UnknownQuantity

class Ions(UnknownQuantity):
    

    def __init__(self, settings):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.ions = list()
        self.r    = None
        self.t    = None


    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, Z0=None, n=None, r=None, t=None, tritium=False):
        """
        Adds a new ion species to the plasma.

        :param str name:        Name by which the ion species will be referred to.
        :param int Z:           Ion charge number.
        :param int iontype:     Method to use for evolving ions in time.
        :param int Z0:          Charge state to populate (used for populating exactly one charge state for the ion).
        :param n:               Ion density (can be either a scalar, 1D array or 2D array, depending on the other input parameters)
        :param numpy.ndarray r: Radial grid on which the input density is defined.
        :param numpy.ndarray t: Time grid on which the input density is defined.
        :param bool tritium:    If ``True``, the ion species is treated as Tritium.
        """
        if (self.r is not None) and (r is not None) and (np.any(self.r != r)):
            raise EquationException("The radial grid must be the same for all ion species.")
        if (self.t is not None) and (t is not None) and (np.any(self.t != t)):
            raise EquationException("The time grid must be the same for all ion species.")

        ion = IonSpecies(settings=self.settings, name=name, Z=Z, ttype=iontype, Z0=Z0, n=n, r=r, t=t, interpr=self.r, interpt=None, tritium=tritium)
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


    def getIon(self, i=-1, name=None):
        """
        Returns the ion species with the specified index.
        """
        if i >= 0: return self.ions[i]
        elif name is not None:
            for i in range(0, len(self.ions)):
                if self.ions[i].getName() == name:
                    return self.ions[i]

            raise EquationException("No ion with name '{}' has been defined.".format(name))
        else:
            raise EquationException("Invalid call to 'getIon()'.")


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
    

    def fromdict(self, data):
        """
        Load settings from the specified dictionary.
        """
        names        = data['names'].split(';')[:-1]
        Z            = data['Z']
        types        = data['types']

        if 'tritiumnames' in data:
            tritiumnames = data['tritiumnames'].split(';')[:-1]
        else:
            tritiumnames = []

        initial    = None
        prescribed = None

        if 'initial' in data:
            initial = data['initial']
        if 'prescribed' in data:
            prescribed = data['prescribed']

        iidx, pidx = 0, 0
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

            tritium = (names[i] in tritiumnames)
            self.addIon(name=names[i], Z=Z[i], iontype=types[i], n=n, r=r, t=t, tritium=tritium)

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this Ions object.
        """
        Z       = self.getCharges()
        itypes  = self.getTypes()
        initial = None
        prescribed = None
        names   = ""
        tritiumnames = ""

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


        data = {
            'names': names,
            'Z': Z,
            'types': itypes
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


    def getFreeElectronDensity(self, t=0):
        n_free = np.zeros( self.r.shape )

        for ion in self.ions:
            for Z0 in range(1,ion.Z + 1):
                if len( ion.n.shape ) == 3:
                    n_free = n_free + Z0 * ion.n[Z0,t,:]
                elif len( ion.n.shape ) == 2:
                    n_free = n_free + Z0 * ion.n[Z0,:]
                
        return n_free, self.r


