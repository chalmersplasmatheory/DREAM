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


    ############################
    # SETTERS
    ############################
    def setFilename(self, filename):
        """
        Set the name of the output file.
        """
        self.filename = filename


    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.filename = data['filename']

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this Output object.
        """
        if verify:
            self.verifySettings()

        data = {
            'filename': self.filename
        }

        return data


    def verifySettings(self):
        """
        Verify that the Output settings are consistent.
        """
        if type(self.filename) != str:
            raise DREAMException("The output file name must be string.")


