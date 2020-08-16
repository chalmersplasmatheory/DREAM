#
# Object handling settings for output from DREAM.
# ################################################

import numpy as np


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
        """
        self.filename = filename


    def setTiming(self, stdout=None, file=None):
        """
        Specifies whether to print timing information and/or include
        it in the output file.
        """
        if stdout is not None:
            self.timingstdout = stdout
        if file is not None:
            self.timingfile = file


    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.filename = data['filename']
        self.timingstdout = data['timingstdout']
        self.timingfile = data['timingfile']

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this Output object.
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


