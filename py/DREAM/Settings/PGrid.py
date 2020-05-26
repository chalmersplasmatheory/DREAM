# Settings for a p (momentum) grid

import numpy as np
from DREAM.DREAMException import DREAMException


class PGrid:

    TYPE_UNIFORM = 1
    
    def __init__(self, name, ttype=1, np=100, pmax=None, data=None):
        """
        Constructor.

          name:  Name of grid (e.g. 'hottailgrid' or 'runawaygrid')
        AND
          ttype: Grid type.
          np:    Number of p grid points.
          pmax:  Maximum value of p.
        OR
          data:  Dictionary containing all of the above settings
                 (except 'ttype' should be called 'pgrid')
        """
        self.name = name

        if data is not None:
            self.fromdict(data)
        else:
            self.setType(ttype=ttype)
            self.setNp(np)

            if pmax is not None:
                self.setPmax(pmax)
            else:
                self.pmax = pmax


    ####################
    # GETTERS
    ####################
    def getNp(self): return self.np
    def getPmax(self): return self.pmax
    def getType(self): return self.type


    ####################
    # SETTERS
    ####################
    def setNp(self, np):
        if np <= 0:
            raise DREAMException("PGrid {}: Invalid value assigned to 'np': {}. Must be > 0.".format(self.name, np))
        elif np == 1:
            print("WARNING: PGrid {}: np = 1. Consider disabling the hot-tail grid altogether.".format(self.name))

        self.np = int(np)


    def setPmax(self, pmax):
        if pmax <= 0:
            raise DREAMException("PGrid {}: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.name, pmax))

        self.pmax = float(pmax)


    def setType(self, ttype):
        """
        Set the type of p grid generator.
        """
        if ttype == self.TYPE_UNIFORM:
            self.type = ttype
        else:
            raise DREAMException("PGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))


    def fromdict(self, data):
        """
        Load this p-grid from the specified dictionary.
        """
        self.type = data['pgrid']
        self.np   = data['np']
        self.pmax = data['pmax']

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this PGrid object.
        """
        if verify:
            self.verifySettings()

        return {
            'pgrid': self.type,
            'np': self.np,
            'pmax': self.pmax
        }


    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type == self.TYPE_UNIFORM:
            if self.np is None or self.np <= 0:
                raise DREAMException("PGrid {}: Invalid value assigned to 'np': {}. Must be > 0.".format(self.name, self.np))
            elif self.pmax is None or self.pmax <= 0:
                raise DREAMException("PGrid {}: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.name, self.pmax))
        else:
            raise DREAMException("PGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))


