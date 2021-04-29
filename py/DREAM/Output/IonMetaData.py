# Object for handling ion meta data (atomic charge and names)

import numpy as np


class IonMetaData:
    

    def __init__(self, data):
        """
        Constructor.

        data: Dictionary containing atomic charge 'Z' and ion names 'names'
              (either as a Python list or DREAM monolithic string list)
        """

        self.Z = [int(Z) for Z in data['Z'][:]]
        self.names = data['names'][:]

        # DREAM monolithic string list? (i.e. a series of strings
        # glued together with ';' in between)
        if ';' in self.names:
            self.names = self.names.split(';')[:-1]


    def __getitem__(self, i):
        """
        IF i IS AN INTEGER
          Returns a tuple consisting of the ion name and
          charge for the specified ion.

        IF i IS A STRING
          Returns the ion charge for the named ion.
        """
        if type(i) == str:
            return self.Z[i]
        else:
            return (self.names[i], self.Z[i])


    def getName(self, i):
        """
        Get the name of the ion with the specified index.
        """
        return self.names[i]


    def getNames(self):
        """
        Get a list of all ion names.
        """
        return self.names


    def getCharge(self, i):
        """
        Get the ion charge for the ion with the specified index.
        """
        return self.Z[i]


    def getCharges(self):
        """
        Get a list of all ion charges.
        """
        return self.Z


