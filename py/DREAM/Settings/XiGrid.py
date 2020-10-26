# Settings for a p (momentum) grid

import numpy as np
from DREAM.DREAMException import DREAMException


TYPE_UNIFORM = 1
TYPE_BIUNIFORM = 2
TYPE_UNIFORM_THETA = 3
TYPE_BIUNIFORM_THETA = 4
    

class XiGrid:
    TYPE_UNIFORM = 1
    TYPE_BIUNIFORM = 2

    def __init__(self, name, ttype=TYPE_UNIFORM, nxi=25, data=None):
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
            self.nxisep = None
            self.xisep  = None
            self.nthetasep = None
            self.thetasep  = None
            


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


    def setBiuniform(self, xisep=None, nxisep = None, nxisep_frac = None,thetasep = None, nthetasep =None, nthetasep_frac=None ):
       
        if xisep is not None:
            self.type = TYPE_BIUNIFORM
            self.xisep = float(xisep)
            if nxisep is not None:
                self.nxisep = int(nxisep)
            elif nxisep_frac is not None:
       	        self.nxisep = int(round(self.nxi * nxisep_frac))
            else:
                raise DREAMException("XiGrid biuniform {}: nxisep or nxisep_frac must be set.")
        elif thetasep is not None:
            self.type = TYPE_BIUNIFORM_THETA
            self.thetasep = float(thetasep)
            if nthetasep is not None:
                self.nthetasep = int(nthetasep)
            elif nthetasep_frac is not None:
                self.nthetasep = int(round(self.nxi*nthetasep_frac))
            else:
                raise DREAMException("XiGrid biuniform theta {}: nthetasep or nthetasep_frac must be set.")
        else:	
            raise DREAMException("XiGrid biuniform  {}: thetasep or xisep must be set.")


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
            self.nxisep = data['nxisep']
            self.xisep  = data['xisep']
        
        elif self.type == TYPE_BIUNIFORM_THETA:
            self.nthetasep = data['nxisep']
            self.thetasep  = data['xisep']
            
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
            data['nxisep'] = self.nxisep
            data['xisep'] = self.xisep
            
        elif self.type == TYPE_BIUNIFORM_THETA:
            data['nxisep'] = self.nthetasep
            data['xisep'] = self.thetasep

        return data

    
    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type in [TYPE_UNIFORM,TYPE_BIUNIFORM,TYPE_UNIFORM_THETA,TYPE_BIUNIFORM_THETA]:
            if self.nxi is None or self.nxi <= 0:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nxi': {}. Must be > 0.".format(self.name, self.nxi))
        else:
            raise DREAMException("XiGrid {}: Unrecognized grid type specified: {}.".format(self.name, self.type))
        if self.type == TYPE_BIUNIFORM:
            if self.nxisep is None or self.nxisep <= 0 or self.nxisep >= self.nxi:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nxisep': {}. Must be > 0 and < nxi.".format(self.name, self.nxisep))
            elif self.xisep is None or self.xisep <= -1 or self.xisep >= 1:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'xisep': {}. Must be > -1 and < 1.".format(self.name, self.xisep))
                
        elif self.type == TYPE_BIUNIFORM_THETA:
            if self.nthetasep is None or self.nthetasep <= 0 or self.nthetasep >= self.nxi:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'nthetasep': {}. Must be > 0 and < nxi.".format(self.name, self.nthetasep))
            elif self.thetasep is None or self.thetasep <= -1 or self.thetasep >= 1:
                raise DREAMException("XiGrid {}: Invalid value assigned to 'thetasep': {}. Must be > -1 and < 1.".format(self.name, self.thetasep))


