
import numpy as np


class ColdElectrons:
    
    TYPE_PRESCRIBED = 1

    def __init__(self):
        pass


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }

        if self.type == self.TYPE_PRESCRIBED:
            data['data'] = {
                'n': self.density,
                'r': self.radialgrid,
                't': self.times
            }
        else:
            raise EquationException("Unrecognized cold electron density type: {}".format(self.type))

        return data
