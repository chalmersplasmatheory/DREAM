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
        'fluid',
        'fluid/conductivity',
        'fluid/Eceff',
        'fluid/GammaAva',
        'fluid/gammaCompton',
        'fluid/gammaDreicer',
        'fluid/lnLambdaC', 'fluid/lnLambdaT',
        'fluid/pCrit',
        'fluid/radiation',
        'fluid/runawayRate',
        'hottail/lnLambda_ee_f1', 'hottail/lnLambda_ee_f2',
        'hottail/lnLambda_ei_f1', 'hottail/lnLambda_ei_f2',
        'hottail/nu_D_f1', 'hottail/nu_D_f2',
        'hottail/nu_s_f1', 'hottail/nu_s_f2',
        'lnLambda',
        'nu_s',
        'nu_D',
        'runaway/lnLambda_ee_f1', 'runaway/lnLambda_ee_f2',
        'runaway/lnLambda_ei_f1', 'runaway/lnLambda_ei_f2',
        'runaway/nu_D_f1', 'runaway/nu_D_f2',
        'runaway/nu_s_f1', 'runaway/nu_s_f2',
        'scalar/radialloss_n_re',
        'transport'
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


    def fromdict(self, data):
        """
        Load these settings from the given dictionary.
        """
        inc = []
        if 'include' in data:
            inc = data['include'].split(';')

        if len(inc) > 0 and inc[-1] == '':
            inc = inc[:-1]

        self.include(inc)


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


