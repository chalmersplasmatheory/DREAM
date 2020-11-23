# Settings for the equation system preconditioner

import numpy as np
from .. DREAMException import DREAMException
from . import EquationSystem


class Preconditioner:
    

    def __init__(self):
        """
        Constructor.
        """
        self.enabled = True
        self.overrides = []


    def fromdict(self, data):
        """
        Load preconditioner settings from a dictionary.
        """
        if 'enabled' in data:
            self.enabled = bool(data['enabled'])

        if 'names' in data:
            if 'equation_scales' not in data:
                raise DREAMException("'names' setting present, but no 'equation_scales' setting found.")
            if 'unknown_scales' not in data:
                raise DREAMException("'names' setting present, but no 'unknown_scales' setting found.")

            names = data['names'].split(';')[:-1]
            for i in range(len(names)):
                escal = data['equation_scales'][i]
                uscal = data['unknown_scales'][i]

                l = {'name': names[i], 'equation_scale': float(escal), 'unknown_scale': float(uscal)}
                overrides.append(l)

            self.overrides = overrides


    def getIndex(self, unknown):
        """
        Returns the index into the 'overrides' list for the
        given unknown. If the returned value is '-1', no
        override exists for the quantity and the default
        scsales are used instead.
        """
        for i in range(0, len(self.overrides)):
            if self.overrides[i]['name'] == unknown:
                return i
        
        return -1


    def set(self, unknown, scale, equation_scale=None):
        """
        Set the scales for one or more unknown quantities.

        :param unknown:        A string or list of strings specifying the name(s) of the quantity/ies to set the scales for.
        :param scale:          Scale to normalize the unknown quantity to.
        :param equation_scale: Scale to normalize the equation of the unknown to. If 'None', it is taken to be the same as 'scale'.
        """
        if equation_scale is None:
            equation_scale = scale

        if type(unknown) == str:
            t = self.getIndex(unknown)
            l = {'name': unknown, 'equation_scale': float(equation_scale), 'unknown_scale': float(scale)}

            if t < 0:
                self.overrides.append(l)
            else:
                self.overrides[t] = l
        elif type(unknown) == list:
            for u in unknown:
                self.set(u, scale=scale, equation_scale=equation_scale)
        else:
            raise DREAMException("Preconditioner.set(): Unrecognized type of parameter 'unknown': {}.".format(type(unknown)))


    def setEnabled(self, enabled=True):
        """
        Enable/disable the physics-based preconditioner.

        :param bool enabled: Indicates whether to enable/disable the preconditioner.
        """
        self.enabled = enabled


    def todict(self):
        """
        Convert this object to a dict.
        """
        data = {'enabled': self.enabled}

        if len(self.overrides) > 0:
            data['names'] = ''
            data['equation_scales'] = []
            data['unknown_scales'] = []

            for u in self.overrides:
                data['names'] += '{};'.format(u['name'])
                data['equation_scales'].append(u['equation_scale'])
                data['unknown_scales'].append(u['unknown_scale'])

        return data


    def verifySettings(self):
        """
        Verify that these settings are consistent.
        """
        if type(self.enabled) is not bool:
            raise DREAMException("Invalid type of option 'enabled': {}. Expected bool.".format(type(self.enabled)))

        for i in range(len(self.overrides)):
            u = self.overrides[i]
            if ('name' not in u) or ('equation_scales' not in u) or ('unknown_scales' not in u):
                raise DREAMException("Incomplete override at index {}.".format(i))
            elif type(u['equation_scale']) is not float:
                raise DREAMException("Invalid type of equation scale for unknown '{}': {}. Expected float.".format(u['name'], type(u['equation_scale'])))
            elif type(u['unknown_scale']) is not float:
                raise DREAMException("Invalid type of scale for unknown '{}': {}. Expected float.".format(u['name'], type(u['unknown_scale'])))

            if u['name'] not in EquationSystem.UNKNOWNS:
                print("WARNING: No unknown quantity with name '{}' is available in DREAM. Preconditioner scale setting will be ignored...".format(u['name']))


