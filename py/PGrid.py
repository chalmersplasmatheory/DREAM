# Settings for a p (momentum) grid

import numpy as np


class PGrid:

    TYPE_UNIFORM = 1
    
    def __init__(self, ttype=1, np=100, pmax=None):
        """
        Constructor.
        """
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
            raise DREAMException("PGrid: Invalid value assigned to 'np': {}. Must be > 0.".format(np))
        elif np == 1:
            print("WARNING: PGrid: np = 1. Consider disabling the hot-tail grid altogether.")

        self.np = int(np)


    def setPmax(self, pmax):
        if pmax <= 0:
            raise DREAMException("PGrid: Invalid value assigned to 'pmax': {}. Must be > 0.".format(pmax))

        self.pmax = float(pmax)


    def setType(self, ttype):
        """
        Set the type of p grid generator.
        """
        if ttype == self.TYPE_UNIFORM:
            self.type = ttype
        else:
            raise DREAMException("PGrid: Unrecognized grid type specified: {}.".format(self.type))


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this PGrid object.
        """
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
                raise DREAMException("PGrid: Invalid value assigned to 'np': {}. Must be > 0.".format(self.np))
            elif self.pmax is None or self.pmax <= 0:
                raise DREAMException("PGrid: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.pmax))
        else:
            raise DREAMException("PGrid: Unrecognized grid type specified: {}.".format(self.type))


