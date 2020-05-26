# Settings for a p (momentum) grid

import numpy as np
from DREAM.DREAMException import DREAMException


TYPE_UNIFORM = 1
    

class XiGrid:

    def __init__(self, name, ttype=1, nxi=25, data=None):
        """
        Constructor.

          name:  Name of grid (e.g. 'hottailgrid' or 'runawaygrid')
        AND
          ttype: Grid type.
          np:    Number of p grid points.
          pmax:  Maximum value of p.
        OR
          data:  Dictionary containing all of the above settings
                 (except 'ttype' should be called 'xigrid')
        """
        self.name = name

        if data is not None:
            self.fromdict(data)
        else:
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
            raise DREAMException("XiGrid {}: Invalid value assigned to 'nxi': {}. Must be > 0.".format(self.name, nxi))

        self.nxi = int(nxi)


    def setType(self, ttype):
        """
        Set the type of xi grid generator.
        """
        if ttype == TYPE_UNIFORM:
            self.type = ttype
        else:
            raise DREAMException("XiGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))


    def fromdict(self, data):
        """
        Load this xi-grid from the specified dictionary.
        """
        self.type = data['xigrid']
        self.nxi  = data['nxi']

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this XiGrid object.
        """
        if verify:
            self.verifySettings()

        return {
            'xigrid': self.type,
            'nxi': self.nxi
        }

    
    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type == TYPE_UNIFORM:
            if self.nxi is None or self.nxi <= 0:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nxi': {}. Must be > 0.".format(self.name, self.np))
        else:
            raise DREAMException("XiGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))


