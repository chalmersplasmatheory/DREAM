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
        'fluid/gammaCompton', 'fluid/gammaDreicer', 'fluid/gammaTritium', 'fluid/gammaHottail',
        'fluid/Lambda_hypres',
        'fluid/lnLambdaC', 'fluid/lnLambdaT',
        'fluid/nusnuDatPstar',
        'fluid/pCrit', 'fluid/pCritHottail',
        'fluid/pStar',
        'fluid/qR0',
        'fluid/radiation',
        'fluid/runawayRate',
        'fluid/Tcold_ohmic',
        'fluid/Tcold_fhot_coll',
        'fluid/Tcold_fre_coll',
        'fluid/Tcold_transport',
        'fluid/Tcold_radiation',
        'fluid/Tcold_binding_energy',
#        'fluid/Tcold_radiationFromNuS',
        'fluid/Tcold_ion_coll',
        'fluid/tIoniz',
        'fluid/W_hot',
        'fluid/W_re',
        'fluid/reioniz_rate', 'fluid/reioniz_vsigma',
        'fluid/Zeff',
        'energy',
        'hottail/Ar', 'hottail/Ap1', 'hottail/Ap2',
        'hottail/Drr', 'hottail/Dpp', 'hottail/Dpx', 'hottail/Dxp', 'hottail/Dxx',
        'hottail/timevaryingb_Ap2',
        'hottail/kinioniz_rate', 'hottail/kinioniz_vsigma',
        'hottail/lnLambda_ee_f1', 'hottail/lnLambda_ee_f2',
        'hottail/lnLambda_ei_f1', 'hottail/lnLambda_ei_f2',
        'hottail/nu_D_f1', 'hottail/nu_D_f2',
        'hottail/nu_s_f1', 'hottail/nu_s_f2',
        'hottail/nu_par_f1', 'hottail/nu_par_f2',
        'hottail/S_ava', 'hottail/S_compton', 'hottail/S_tritium',
        'hottail/synchrotron_f1', 'hottail/synchrotron_f2', 'hottail/bremsstrahlung_f1', 
        'lnLambda',
        'nu_s',
        'nu_D',
        'runaway/Ar', 'runaway/Ap1', 'runaway/Ap2',
        'runaway/Drr', 'runaway/Dpp', 'runaway/Dpx', 'runaway/Dxp', 'runaway/Dxx',
        'runaway/timevaryingb_Ap2',
        'runaway/kinioniz_rate', 'runaway/kinioniz_vsigma',
        'runaway/lnLambda_ee_f1', 'runaway/lnLambda_ee_f2',
        'runaway/lnLambda_ei_f1', 'runaway/lnLambda_ei_f2',
        'runaway/nu_D_f1', 'runaway/nu_D_f2',
        'runaway/nu_s_f1', 'runaway/nu_s_f2',
        'runaway/nu_par_f1', 'runaway/nu_par_f2',
        'runaway/S_ava', 'runaway/S_compton', 
        'runaway/synchrotron_f1', 'runaway/synchrotron_f2', 'runaway/bremsstrahlung_f1', 
        'scalar',
        'scalar/ablationDrift',
        'scalar/E_mag',
        'scalar/L_i',
        'scalar/L_i_flux',
        'scalar/l_i',
        'scalar/radialloss_n_re',
        'scalar/energyloss_T_cold',
        'scalar/radialloss_f_re',
        'scalar/radialloss_f_hot',
        'scalar/energyloss_f_re',
        'scalar/energyloss_f_hot',
        'scalar/tIoniz',
        'scalar/Ypdot',
        'ripple',
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
