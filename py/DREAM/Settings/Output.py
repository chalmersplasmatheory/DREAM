#
# Object handling settings for output from DREAM.
# ################################################

import numpy as np

from .. DREAMException import DREAMException


class Output:
    
    
    def __init__(self, filename='output.h5'):
        """
        Constructor.
        """
        self.filename = filename
        self.timingstdout = False
        self.timingfile = True


    ############################
    # SETTERS
    ############################
    def setFilename(self, filename):
        """
        Set the name of the output file.

        :param str filename: Name of output file to store simulation data to.
        """
        self.filename = filename


    def setTiming(self, stdout=None, file=None):
        """
        Specifies whether to print timing information and/or include
        it in the output file.

        :param bool stdout: If ``True``, prints details about execution time to ``stdout`` at the end of the simulation.
        :param bool file:   If ``True``, stores details about execution time in the output file after the simulation.
        """
        if stdout is not None:
            self.timingstdout = stdout
        if file is not None:
            self.timingfile = file


    def fromdict(self, data):
        """
        Load settings from the given dictionary.

        :param dict data: Dictionary to load settings from.
        """
        self.filename = data['filename']
        self.timingstdout = bool(data['timingstdout'])
        self.timingfile = bool(data['timingfile'])

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this Output object.

        :param bool verify: If ``True``, verifies the settings of this object before generating the dictionary.
        """
        if verify:
            self.verifySettings()

        data = {
            'filename': self.filename,
            'timingfile': self.timingfile,
            'timingstdout': self.timingstdout
        }

        return data


    def verifySettings(self):
        """
        Verify that the Output settings are consistent.
        """
        if type(self.filename) != str:
            raise DREAMException("The output file name must be string.")
        elif type(self.timingfile) != bool:
            raise DREAMException("The option 'timingfile' must be a bool.")
        elif type(self.timingstdout) != bool:
            raise DREAMException("The option 'timingstdout' must be a bool.")


