# Settings for a p (momentum) grid

import numpy as np
from DREAM.DREAMException import DREAMException
from DREAM.Settings.Equations.EquationException import EquationException


TYPE_UNIFORM = 1
TYPE_BIUNIFORM = 2
TYPE_UNIFORM_THETA = 3
TYPE_BIUNIFORM_THETA = 4
TYPE_CUSTOM = 5


class XiGrid:
    def __init__(self, name, ttype=TYPE_UNIFORM, nxi=0, data=None):
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

        self.nxisep = None
        self.nxisep_frac = None
        self.xisep  = None
        self.nthetasep = None
        self.nthetasep_frac = None
        self.thetasep  = None

        self.xi_f = None

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
        self.nxi = int(nxi)


    def setBiuniform(self, xisep=None, nxisep = None, nxisep_frac = None,thetasep = None, nthetasep =None, nthetasep_frac=None ):
       
        if xisep is not None:
            self.type = TYPE_BIUNIFORM
            self.xisep = float(xisep)
            if nxisep is not None:
                self.nxisep = int(nxisep)
                self.nxisep_frac = None
            elif nxisep_frac is not None:
                self.nxisep = None
                self.nxisep_frac = nxisep_frac
            else:
                raise DREAMException("XiGrid biuniform {}: nxisep or nxisep_frac must be set.")
        elif thetasep is not None:
            self.type = TYPE_BIUNIFORM_THETA
            self.thetasep = float(thetasep)
            if nthetasep is not None:
                self.nthetasep = int(nthetasep)
                self.nthetasep_frac = None
            elif nthetasep_frac is not None:
                self.nthetasep = None
                self.nthetasep_frac = nthetasep_frac
            else:
                raise DREAMException("XiGrid biuniform theta {}: nthetasep or nthetasep_frac must be set.")
        else:	
            raise DREAMException("XiGrid biuniform  {}: thetasep or xisep must be set.")

    def setCustomGridPoints(self, xi_f):
        """
        Set an arbitrary custom grid point distribution
        on the pitch flux grid (i.e. the locations of
        the cell edges). This overrides the grid resolution
        'nxi', which will be taken as the number of cells
        described by the prescribed grid points.

        :param float xi_f: List of pitch flux grid points
        """
        self.type = TYPE_CUSTOM
        if type(xi_f)==list:
            xi_f = np.array(xi_f)
        if np.size(xi_f)<2:
            raise EquationException("XiGrid: Custom grid point vector 'xi_f' must have size 2 or greater.")
        for i in range(np.size(xi_f)-1):
            if xi_f[i+1]<xi_f[i]:
                raise EquationException("XiGrid: Custom grid points 'xi_f' must be an array of increasing numbers.")
        if np.min(xi_f)!=-1 or np.max(xi_f)!=1:
            raise EquationException("XiGrid: Custom pitch grid must span [-1,1].")
        self.xi_f = xi_f
        if self.nxi != 0:
            print("*WARNING* XiGrid: Prescibing custom pitch grid overrides 'nxi'.")
        self.nxi = np.size(self.xi_f) - 1


    def setType(self, ttype):
        """
        Set the type of xi grid generator.
        """
        if ttype in [TYPE_UNIFORM,TYPE_BIUNIFORM,TYPE_UNIFORM_THETA,TYPE_BIUNIFORM_THETA]:
            self.type = ttype
        else:
            raise DREAMException("XiGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))


    def fromdict(self, data):
        """
        Load this xi-grid from the specified dictionary.
        """
        self.type = data['xigrid']
        self.nxi  = data['nxi']
        if self.type == TYPE_BIUNIFORM:
            self.xisep  = data['xisep']
            if 'nxisep' in data:
                self.nxisep = data['nxisep']
            elif 'nxisep_frac' in data:
                self.nxisep_frac = data['nxisep_frac']
        
        elif self.type == TYPE_BIUNIFORM_THETA:
            self.thetasep  = data['xisep']
            if 'nxisep' in data:
                self.nthetasep = data['nxisep']
            elif 'nxisep_frac' in data:
                self.nthetasep_frac = data['nxisep_frac']
        elif self.type == TYPE_CUSTOM:
            self.xi_f = data['xi_f']
            
        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this XiGrid object.
        """
        if verify:
            self.verifySettings()

        data = { 
            'xigrid': self.type, 
            'nxi': self.nxi,
        }
        if self.type == TYPE_BIUNIFORM:
            if self.nxisep is not None:
                data['nxisep'] = self.nxisep
            elif self.nxisep_frac is not None:
                #data['nxisep'] = int(round(self.nxi*self.nxisep_frac))
                data['nxisep_frac'] = self.nxisep_frac
            else:
                raise DREAMException("XiGrid {}: Neither 'nxisep' nor 'nxisep_frac' have been specified.".format(self.name))

            data['xisep'] = self.xisep
            
        elif self.type == TYPE_BIUNIFORM_THETA:
            if self.nthetasep is not None:
                data['nxisep'] = self.nthetasep
            elif self.nthetasep_frac is not None:
                #data['nxisep'] = int(round(self.nxi*self.nthetasep_frac))
                data['nxisep_frac'] = self.nthetasep_frac

            data['xisep'] = self.thetasep
        elif self.type == TYPE_CUSTOM:
            data['xi_f'] = self.xi_f

        return data

    
    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type in [TYPE_UNIFORM,TYPE_BIUNIFORM,TYPE_UNIFORM_THETA,TYPE_BIUNIFORM_THETA,TYPE_CUSTOM]:
            if self.nxi is None or self.nxi <= 0:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nxi': {}. Must be > 0.".format(self.name, self.nxi))
        else:
            raise DREAMException("XiGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))
        if self.type == TYPE_BIUNIFORM:
            if self.nxisep is not None and (self.nxisep <= 0 or self.nxisep >= self.nxi):
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nxisep': {}. Must be > 0 and < nxi.".format(self.name, self.nxisep))
            if self.nxisep_frac is not None and (self.nxisep_frac <= 0 or self.nxisep_frac >= 1):
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nxisep_frac': {}. Must be > 0 and < 1.".format(self.name, self.nxisep))
            elif self.nxisep is None and self.nxisep_frac is None:
                raise DREAMException("XiGrid {}: Neither 'nxisep' nor 'nxisep_frac' have been specified.".format(self.name))
            elif self.xisep is None or self.xisep <= -1 or self.xisep >= 1:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'xisep': {}. Must be > -1 and < 1.".format(self.name, self.xisep))
                
        elif self.type == TYPE_BIUNIFORM_THETA:
            if self.nthetasep is not None and (self.nthetasep <= 0 or self.nthetasep >= self.nxi):
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nthetasep': {}. Must be > 0 and < nxi.".format(self.name, self.nthetasep))
            elif self.nthetasep_frac is not None and (self.nthetasep_frac <= 0 or self.nthetasep_frac >= 1):
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nthetasep_frac': {}. Must be > 0 and < 1.".format(self.name, self.nthetasep_frac))
            elif self.nthetasep is None and self.nthetasep_frac is None:
                raise DREAMException("XiGrid {}: Neither 'nthetasep' nor 'nthetasep_frac' have been specified.".format(self.name))
            elif self.thetasep is None or self.thetasep <= 0 or self.thetasep >= np.pi:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'thetasep': {}. Must be > 0 and < pi.".format(self.name, self.thetasep))


