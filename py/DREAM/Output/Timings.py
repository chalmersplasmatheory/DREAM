# Class for handling timing output


import numpy as np
from ..DREAMException import DREAMException
from ..DataObject import DataObject

class Timings:
    

    def __init__(self, timings=None, output=None):
        """
        Constructor.

        timings: Dictionary containing timings information from the DREAM kernel.
        output:  Parent DREAMOutput object.
        """
        self.timings = {}
        self.output  = output
        self.descriptions = {}
        self.subtimers = []

        if timings is not None:
            self.loadTimingInformation(timings)


    def __contains__(self, item):
        """
        Overrides the Python 'in' operator.
        """
        return (item in self.__dict__)


    def __getitem__(self, index):
        """
        Direct access by name to the timing information.
        """
        return self.__dict__[index]


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()

    
    def __str__(self):
        """
        Convert this object to a string.
        """
        s = ""
        for key in self.timings:
            unit, t = self.formatTimeAndUnit(self.timings[key])
            s += '  {:30s}  -- {} {}\n'.format(self.descriptions[key], t, unit)

        for key in self.subtimers:
            unit, t = self.formatTimeAndUnit(self[key].getTotal())
            s += '  {:30s}  -- {} {}\n'.format(key, t, unit)

        return s


    def formatTimeAndUnit(self, time):
        """
        Format the given time and provide an appropriate unit.
        """
        # Convert to ms
        unit = 'ms'
        t = time / 1000

        if t > 1000:
            t /= 1000
            unit = 's'

        return unit, t


    def loadTimingInformation(self, timings):
        """
        Load timing information from the given dict.
        """
        tim = {}

        for key in timings:
            if key[-2:] == '@@':
                self.descriptions[key[:-2]] = timings[key]['desc']
            elif type(timings[key]) == float:
                setattr(self, key, timings[key])
                tim[key] = timings[key]
            elif type(timings[key]) == np.ndarray:
                setattr(self, key, timings[key][0])
                tim[key] = timings[key][0]
            elif type(timings[key]) == DataObject:
                setattr(self, key, timings[key][:][0])
                tim[key] = timings[key][:][0]
            elif type(timings[key]) == dict:
                setattr(self, key, Timings(timings[key], output=self.output))
                self.subtimers.append(key)
            else:
                raise DREAMException("Unrecognized type of member '{}' of timings information: {}.".format(key, type(timings[key])))

        self.timings   = tim
        self.subtimers = sorted(self.subtimers)
        

    def getTotal(self):
        """
        Get total simulation time.
        """
        if 'total' in self.timings:
            return self.timings['total']
        else:
            t = 0
            for key in self.timings:
                if type(self.timings[key]) == float:
                    t += self.timings[key]
            
            return t
                

