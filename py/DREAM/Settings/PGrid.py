# Settings for a p (momentum) grid

import numpy as np
from DREAM.DREAMException import DREAMException


class PGrid:

    TYPE_UNIFORM = 1
    TYPE_BIUNIFORM = 2
    
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

        self.npsep = None
        self.npsep_frac = None
        self.psep  = None

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

    def setBiuniform(self, psep, npsep = None, npsep_frac = None):
        self.type = self.TYPE_BIUNIFORM
        self.psep = psep
        if npsep is not None:
            self.npsep = npsep
            self.npsep_frac = None
        elif npsep_frac is not None:
            self.npsep = None
            self.npsep_frac = npsep_frac
        else:
            raise DREAMException("PGrid biuniform {}: npsep or npsep_frac must be set.")


    def setType(self, ttype):
        """
        Set the type of p grid generator.
        """
        if ttype == self.TYPE_UNIFORM or ttype == self.TYPE_BIUNIFORM:
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

        if self.type == self.TYPE_BIUNIFORM:
            if 'npsep' in data:
                self.npsep = data['npsep']
            if 'npsep_frac' in data:
                self.npsep_frac = data['npsep_frac']

            self.psep  = data['psep']
        
        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this PGrid object.
        """
        if verify:
            self.verifySettings()


        data = { 
            'pgrid': self.type, 
            'np': self.np,
            'pmax': self.pmax,
        }
        if self.type == self.TYPE_BIUNIFORM:
            if self.npsep is not None:
                data['npsep'] = self.npsep
            elif self.npsep_frac is not None:
                data['npsep_frac'] = self.npsep_frac

            data['psep'] = self.psep

        return data


    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type == self.TYPE_UNIFORM or self.type == self.TYPE_BIUNIFORM:
            if self.np is None or self.np <= 0:
                raise DREAMException("PGrid {}: Invalid value assigned to 'np': {}. Must be > 0.".format(self.name, self.np))
            elif self.pmax is None or self.pmax <= 0:
                raise DREAMException("PGrid {}: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.name, self.pmax))
        else:
            raise DREAMException("PGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))
        if self.type == self.TYPE_BIUNIFORM:
            if self.npsep is not None and (self.npsep <= 0 or self.npsep >= self.np):
                raise DREAMException("PGrid {}: Invalid value assigned to 'npsep': {}. Must be > 0 and < np.".format(self.name, self.npsep))
            elif self.npsep_frac is not None and (self.npsep_frac <= 0 or self.npsep_frac >= 1):
                raise DREAMException("PGrid {}: Invalid value assigned to 'npsep_frac': {}. Must be > 0 and < np.".format(self.name, self.npsep_frac))
            elif self.npsep is None and self.npsep_frac is None:
                raise DREAMException("PGrid {}: Neither 'npsep' nor 'npsep_frac' have been set.".format(self.name))
            elif self.psep is None or self.psep <= 0 or self.psep >= self.pmax:
                raise DREAMException("PGrid {}: Invalid value assigned to 'psep': {}. Must be > 0 and < pmax.".format(self.name, self.psep))
        

