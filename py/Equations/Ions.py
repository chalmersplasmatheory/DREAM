
import matplotlib.pyplot as plt
import numpy as np
from Equations.EquationException import EquationException
from Equations.IonSpecies import IonSpecies


class Ions:
    

    def __init__(self):
        """
        Constructor.
        """
        self.ions = list()
        self.r    = None
        self.t    = None


    def addIon(self, name, Z, iontype=0, n=None, r=None, t=None):
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
        if (self.r is not None) and (r is not None) and (self.r != r):
            raise EquationException("The radial grid must be the same for all ion species.")
        if (self.t is not None) and (t is not None) and (self.t != t):
            raise EquationException("The time grid must be the same for all ion species.")

        ion = IonSpecies(name=name, Z=Z, ttype=iontype, n=n, r=r, t=t, interpr=self.r, interpt=None)
        self.ions.append(ion)

        self.r = ion.getR()
        if ion.getT() is not None:
            self.t = ion.getT()


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
    

    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this Ions object.
        """
        Z       = self.getCharges()
        itypes  = self.getTypes()
        initial = None
        prescribed = None

        for ion in self.ions:
            if ion.getT() is None:
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
        print("WARNING: Not properly verifying ion settings.");
        return True



