# Handler for "other" quantities, such as collision frequencies,
# bounce averages etc.

import numpy as np
from .. DREAMException import DREAMException


class OtherQuantities:
    
    
    # Here, we keep a list of the possible settings found in DREAM.
    # This allows to check the input the user gives, and emit warnings
    # if the user specifies an unrecognized quantity.
    QUANTITIES = [
        'all',
        'nu_s',
        'nu_D',
        'hottail/nu_s_f1', 'hottail/nu_s_f2',
        'runaway/nu_s_f1', 'runaway/nu_s_f2',
        'hottail/nu_D_f1', 'hottail/nu_D_f2',
        'runaway/nu_D_f1', 'runaway/nu_D_f2',
        'fluid/Eceff',
        'fluid/GammaAva',
        'fluid/lnLambdaC', 'fluid/lnLambdaT'
    ]

    def __init__(self):
        """
        Constructor.
        """
        self._include = list()

    
    def include(self, *args):
        """
        Include one or more "other" quantities in the output.
        """
        for a in args:
            if type(a) == list:
                self.include(*a)
            elif type(a) == str:
                if a not in self.QUANTITIES:
                    print("WARNING: Unrecognized other quantity '{}'. Is it perhaps misspelled?".format(a))

                self._include.append(a)
            else:
                raise DREAMException("other: Unrecognized type of argument: '{}'.".format(type(a)))


    def todict(self, verify=True):
        """
        Returns a dict representing the settings in this object.
        """
        if verify:
            self.verifySettings()

        if len(self._include) == 0:
            return {}
        else:
            return {'include': ';'.join(self._include)}


    def verifySettings(self):
        """
        Verify that these settings are consistent.
        """
        pass


