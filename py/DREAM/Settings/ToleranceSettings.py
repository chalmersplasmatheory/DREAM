# Settings for tolerances of computed quantities

import numpy as np
from .. DREAMException import DREAMException
from . import EquationSystem


class ToleranceSettings:
    
    def __init__(self):
        """
        Constructor.
        """
        # General absolute and relative tolerances (applied to
        # unknowns without explicit tolerances)
        self.reltol = 1e-6
        self.overrides = []


    def disable(self, unknown=None):
        """
        Disable tolerance checking for the given unknowns. If 'unknown'
        is ``None``, tolerance checking is disabled for all unknowns.
        This can be reset for each unknown by explicitly calling
        :automethod:`set` for them.
        """
        if unknown is None:
            self.reltol = 0.0
            self.overrides = []
        elif type(unknown) == str:
            self.set(unknown, abstol=0, reltol=0)
        elif type(unknown) == list:
            for u in unknown:
                self.set(u, abstol=0, reltol=0)
        else:
            raise DREAMException("Unrecognized type of parameter 'unknown': {}.".format(type(unknown)))


    def fromdict(self, data):
        """
        Load tolerance settings from a dictionary.
        """
        if 'reltol' in data:
            r = data['reltol']
            if type(r) == float: self.reltol = r
            else: self.reltol = float(r[0])

        overrides = []

        if 'names' in data:
            if 'abstols' not in data:
                raise DREAMException("'names' setting present, but no 'abstols' setting found.")
            if 'reltols' not in data:
                raise DREAMException("'names' setting present, but no 'reltols' setting found.")
                
            names = data['names'].split(';')[:-1]
            for i in range(len(names)):
                atol = data['abstols'][i]
                rtol = data['reltols'][i]

                l = {'name': names[i], 'abstol': float(atol), 'reltol': float(rtol)}
                overrides.append(l)

            self.overrides = overrides

    
    def getIndex(self, unknown):
        """
        Returns the index into the 'overrides' list for the
        given unknown. If the returned value is '-1', no
        override exists for the quantity and the general
        absolute and relative tolerances are used instead.
        """
        for i in range(0, len(self.overrides)):
            if self.overrides[i]['name'] == unknown:
                return i
        
        return -1


    def set(self, unknown=None, abstol=None, reltol=None):
        """
        Set the absolute and relative tolerance of one or more unknown
        quantities.

        :param unknown: A string or list of strings specifying the name(s) of the quantity/ies to set the tolerances for. If 'None', the tolerances are applied to all unknowns.
        :param float abstol: Absolut tolerance to set for unknown.
        :param float reltol: Relative tolerance to set for unknown.
        """
        if unknown is None:
            if reltol is None:
                raise Exception('If no unknown is specified, the default relative tolerance to use for the system must be specified.')

            self.reltol = float(reltol)
            self.overrides = []
        elif type(unknown) == str:
            t = self.getIndex(unknown)

            # Append to list or overwrite existing element?
            if t < 0:
                if abstol is None: abstol = 0
                if reltol is None: reltol = self.reltol

                self.overrides.append({'name': unknown, 'abstol': float(abstol), 'reltol': float(reltol)})
            else:
                if abstol is not None:
                    self.overrides[t] = float(abstol)
                if reltol is not None:
                    self.overrides[t] = float(reltol)
        elif type(unknown) == list:
            for u in unknown:
                self.set(u, abstol=abstol, reltol=reltol)
        else:
            raise DREAMException("ToleranceSettings.set(): Unrecognized type of parameter 'unknown': {}.".format(type(unknown)))


    def todict(self):
        """
        Convert this object to a dict.
        """
        data = {'reltol': self.reltol}

        if len(self.overrides) > 0:
            data['names'] = ''
            data['abstols'] = []
            data['reltols'] = []
            for u in self.overrides:
                data['names'] += '{};'.format(u['name'])
                data['abstols'].append(u['abstol'])
                data['reltols'].append(u['reltol'])

        return data


    def verifySettings(self):
        """
        Verify that these settings are consistent.
        """
        if type(self.reltol) is not float:
            raise DREAMException("Invalid type of general relative tolerance: {}. Expected float.".format(type(self.reltol)))


        for i in range(len(self.overrides)):
            u = self.overrides[i]
            if ('name' not in u) or ('abstol' not in u) or ('reltol' not in u):
                raise DREAMException("Incomplete override at index {}.".format(i))
            elif type(u['abstol']) is not float:
                raise DREAMException("Invalid type of absolute tolerance for unknown '{}': {}. Expected float.".format(u['name'], type(u['abstol'])))
            elif type(u['reltol']) is not float:
                raise DREAMException("Invalid type of relative tolerance for unknown '{}': {}. Expected float.".format(u['name'], type(u['reltol'])))

            if u['name'] not in EquationSystem.UNKNOWNS:
                print("WARNING: No unknown quantity with name '{}' is available in DREAM. Tolerance setting will be ignored...".format(u['name']))


