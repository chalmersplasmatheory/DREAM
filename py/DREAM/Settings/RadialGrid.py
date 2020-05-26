#
# RadialGrid settings object.
########################################

import numpy as np
from DREAM.DREAMException import DREAMException


TYPE_CYLINDRICAL = 1
TYPE_ANALYTIC_TOROIDAL = 2


class RadialGrid:
    

    def __init__(self, ttype=1):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        # Cylindrical settings
        self.a  = 0.0
        self.B0 = 0.0
        self.nr = int(0)
        self.r0 = 0.0


    #######################
    # SETTERS
    #######################
    def setB0(self, B0):
        if B0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'B0'.")
        
        self.B0 = float(B0)


    def setInnerRadius(self, r0):
        if r0 < 0:
            raise DREAMException("RadialGrid: Invalid value assigned to innermost radius 'r0': {}".format(r0))

        self.r0 = r0


    def setMinorRadius(self, a):
        if a <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(a))

        self.a = float(a)


    def setNr(self, nr):
        if nr <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'nr': {}".format(nr))

        self.nr = int(nr)


    def setType(self, ttype):
        if ttype == TYPE_CYLINDRICAL:
            self.type = ttype
        elif ttype == TYPE_ANALYTIC_TOROIDAL:
            #self.type = ttype
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(ttype))

    
    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.type = data['type']

        if self.type == TYPE_CYLINDRICAL:
            self.a = data['a']
            self.B0 = data['B0']
            self.nr = data['nr']
            self.r0 = data['r0']
        elif self.type == TYPE_ANALYTICAL_TOROIDAL:
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))


    def todict(self, verify=True):
        """
        Returns the settings in this object as a Python dictionary.
        """
        if verify:
            self.verifySettings()

        data = {
            'type': self.type
        }

        if self.type == TYPE_CYLINDRICAL:
            data['a'] = self.a
            data['B0'] = self.B0
            data['nr'] = self.nr
            data['r0'] = self.r0
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        return data
        
            
    def verifySettings(self):
        """
        Verfiy that the RadialGrid settings are consistent.
        """
        if self.type == TYPE_CYLINDRICAL:
            if self.a is None or self.a <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(self.a))
            elif self.B0 is None or self.B0 <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to 'B0': {}".format(self.B0))
            elif self.r0 is None or self.r0 < 0:
                raise DREAMException("RadialGrid: Invalid value assigned to innermost simulated radius 'r0': {}".format(self.r0))

            if self.nr <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned 'nr': {}. Must be > 0.".format(self.nr))

            if self.r0 >= self.a:
                raise DREAMException("RadialGrid: 'r0' must be strictly less than 'a'.")
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))
        

