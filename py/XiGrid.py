# Settings for a p (momentum) grid

import numpy as np


class XiGrid:

    TYPE_UNIFORM = 1
    
    def __init__(self, ttype=1, nxi=25):
        """
        Constructor.
        """
        self.setType(ttype=ttype)
        self.setNxi(nxi)


    ####################
    # GETTERS
    ####################
    def getNxi(self): return self.nxi
    def getType(self): return self.type


    ####################
    # SETTERS
    ####################
    def setNxi(self, nxi):
        if nxi <= 0:
            raise DREAMException("XiGrid: Invalid value assigned to 'nxi': {}. Must be > 0.".format(nxi))

        self.nxi = int(nxi)


    def setType(self, ttype):
        """
        Set the type of xi grid generator.
        """
        if ttype == TYPE_UNIFORM:
            self.type = ttype
        else:
            raise DREAMException("XiGrid: Unrecognized grid type specified: {}.".format(self.type))


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this XiGrid object.
        """
        return {
            'xigrid': self.type,
            'nxi': self.nxi
        }

    
    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type == self.TYPE_UNIFORM:
            if self.np is None or self.np <= 0:
                raise DREAMException("XiGrid: Invalid value assigned to 'np': {}. Must be > 0.".format(self.np))
            elif self.pmax is None or self.pmax <= 0:
                raise DREAMException("XiGrid: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.pmax))
        else:
            raise DREAMException("XiGrid: Unrecognized grid type specified: {}.".format(self.type))


