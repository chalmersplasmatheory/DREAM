
import matplotlib.pyplot as plt
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.IonSpecies import IonSpecies, IONS_PRESCRIBED

class Ions:
    

    def __init__(self):
        """
        Constructor.
        """
        self.ions = list()
        self.r    = None
        self.t    = None


    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, n=None, r=None, t=None):
        """
        Adds a new ion species to the plasma.

        name:    Name by which the ion species will be referred to.
        Z:       Ion charge number.
        iontype: Method to use for evolving ions in time.
        n:       Ion density (can be either a scalar, 1D array or
                 2D array, depending on the other input parameters)
        r:       Radial grid on which the input density is defined.
        t:       Time grid on which the input density is defined.
        """
        if (self.r is not None) and (r is not None) and (np.any(self.r != r)):
            raise EquationException("The radial grid must be the same for all ion species.")
        if (self.t is not None) and (t is not None) and (np.any(self.t != t)):
            raise EquationException("The time grid must be the same for all ion species.")

        ion = IonSpecies(name=name, Z=Z, ttype=iontype, n=n, r=r, t=t, interpr=self.r, interpt=None)
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
        if i > 0: return self.ions[i]
        elif name is not None:
            for i in range(0, len(self.ions)):
                if self.ions[i].getName() == name:
                    return self.ions[i]

            raise EquationException("No ion with name '{}' has been defined.".format(name))
        else:
            raise EquationException("Invalid call to 'getIon()'.")


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
        names = data['names'].split(';')[:-1]
        Z     = data['Z']
        types = data['types']

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
                t = initial['t']
                iidx += Z[i]+1

            self.addIon(name=names[i], Z=Z[i], iontype=types[i], n=n, r=r, t=t)

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

        for ion in self.ions:
            names += '{};'.format(ion.getName())

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
