# Settings for a p (momentum) grid

import numpy as np
from DREAM.DREAMException import DREAMException
from DREAM.Settings.Equations.EquationException import EquationException

TYPE_UNIFORM = 1
TYPE_BIUNIFORM = 2
TYPE_CUSTOM = 3

class PGrid:


    def __init__(self, name, ttype=1, np=0, pmax=None, data=None):
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
        self.p_f = None
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
        if np == 1:
            print("WARNING: PGrid {}: np = 1. Consider disabling the hot-tail grid altogether.".format(self.name))
        self.np = int(np)


    def setPmax(self, pmax):
        if pmax <= 0:
            raise DREAMException("PGrid {}: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.name, pmax))

        self.pmax = float(pmax)

    def setBiuniform(self, psep, npsep = None, npsep_frac = None):
        self.type = TYPE_BIUNIFORM
        self.psep = psep
        if npsep is not None:
            self.npsep = npsep
            self.npsep_frac = None
        elif npsep_frac is not None:
            self.npsep = None
            self.npsep_frac = npsep_frac
        else:
            raise DREAMException("PGrid biuniform {}: npsep or npsep_frac must be set.")

    def setCustomGridPoints(self, p_f):
        """
        Set an arbitrary custom grid point distribution
        on the momentum flux grid (i.e. the locations of
        the cell edges). This overrides the grid resolution
        'np' which will be taken as the number of cells
        described by the prescribed grid points, and 'pmax'
        which will be taken as the largest element in p_f

        :param float p_f: List of momentum flux grid points
        """
        self.type = TYPE_CUSTOM
        if type(p_f)==list:
            p_f = np.array(p_f)
        if np.size(p_f)<2:
            raise EquationException("PGrid: Custom grid point vector 'p_f' must have size 2 or greater.")
        for i in range(np.size(p_f)-1):
            if p_f[i+1]<p_f[i]:
                raise EquationException("PGrid: Custom grid points 'p_f' must be an array of increasing numbers.")
        if np.min(p_f)!=0:
            raise EquationException("PGrid: Custom momentum grid must include 0.")
        self.p_f = p_f
        if self.np != 0:
            print("*WARNING* PGrid: Prescibing custom momentum grid overrides 'np'.")
        self.np = np.size(self.p_f) - 1

        if self.pmax is not None:
            print("*WARNING* PGrid: Prescibing custom momentum grid overrides 'pmax'.")
        self.pmax = float(p_f[-1])        



    def setType(self, ttype):
        """
        Set the type of p grid generator.
        """
        if ttype == TYPE_UNIFORM or ttype == TYPE_BIUNIFORM:
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

        if self.type == TYPE_BIUNIFORM:
            if 'npsep' in data:
                self.npsep = data['npsep']
            if 'npsep_frac' in data:
                self.npsep_frac = data['npsep_frac']

            self.psep  = data['psep']
        elif self.type == TYPE_CUSTOM:
            self.p_f = data['p_f']

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
        if self.type == TYPE_BIUNIFORM:
            if self.npsep is not None:
                data['npsep'] = self.npsep
            elif self.npsep_frac is not None:
                data['npsep_frac'] = self.npsep_frac

            data['psep'] = self.psep
        elif self.type == TYPE_CUSTOM:
            data['p_f'] = self.p_f
        return data


    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type == TYPE_UNIFORM or self.type == TYPE_BIUNIFORM or self.type == TYPE_CUSTOM:
            if self.np is None or self.np <= 0:
                raise DREAMException("PGrid {}: Invalid value assigned to 'np': {}. Must be > 0.".format(self.name, self.np))
            elif self.pmax is None or self.pmax <= 0:
                raise DREAMException("PGrid {}: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.name, self.pmax))
        else:
            raise DREAMException("PGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))
        if self.type == TYPE_BIUNIFORM:
            if self.npsep is not None and (self.npsep <= 0 or self.npsep >= self.np):
                raise DREAMException("PGrid {}: Invalid value assigned to 'npsep': {}. Must be > 0 and < np.".format(self.name, self.npsep))
            elif self.npsep_frac is not None and (self.npsep_frac <= 0 or self.npsep_frac >= 1):
                raise DREAMException("PGrid {}: Invalid value assigned to 'npsep_frac': {}. Must be > 0 and < np.".format(self.name, self.npsep_frac))
            elif self.npsep is None and self.npsep_frac is None:
                raise DREAMException("PGrid {}: Neither 'npsep' nor 'npsep_frac' have been set.".format(self.name))
            elif self.psep is None or self.psep <= 0 or self.psep >= self.pmax:
                raise DREAMException("PGrid {}: Invalid value assigned to 'psep': {}. Must be > 0 and < pmax.".format(self.name, self.psep))
        

