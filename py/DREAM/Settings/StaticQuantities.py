# Handler for "static" quantities

from .. DREAMException import DREAMException


class StaticQuantities:
    

    QUANTITIES = [
        "boltzmann/deltaAngleKernel_hot_hot",
        "boltzmann/deltaAngleKernel_hot_re",
        "boltzmann/deltaAngleKernel_re_re"
    ]


    def __init__(self):
        """
        Constructor.
        """
        self._include = list()


    def include(self, *args):
        """
        Include one or more "static" quantities in the output.
        """
        for a in args:
            if type(a) == list:
                self.include(*a)
            elif type(a) == str:
                if a not in self.QUANTITIES:
                    print(f"WARNING: Unrecognized static quantity '{a}'. Is it perhaps misspelled?")

                self._include.append(a)
            else:
                raise DREAMException(f"static: Unrecognized type of argument: '{type(a)}'.")


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
        Returns the dict representing the settings in this object.
        """
        if verify:
            self.verifySettings()

        if len(self._include) == 0:
            return {}
        else:
            return {'include': ';'.join(self._include) + ';'}


    def verifySettings(self):
        """
        Verify that these settings are consistent.
        """
        pass


