# Class for handling code information output

from ..DataObject import DataObject


class Code:
    

    def __init__(self, code=None, output=None):
        """
        Constructor.
        """
        self.commit = None
        self.refspec = None
        self.datetime_commit = None
        self.datetime_simulation = None

        if code is not None:
            self.loadCodeInfo(code)


    def loadCodeInfo(self, code):
        """
        Load code information.
        """
        self.commit = code['commit']
        self.refspec = code['refspec']
        self.datetime_commit = code['datetime_commit']
        self.datetime_simulation = code['datetime_simulation']
    

    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        s  = "DREAM output generated with version\n"
        s +=f"  {self.refspec}\n"
        s +=f"  SHA1 {self.commit}\n"
        s +=f"  commit made on {self.datetime_commit}\n"

        return s


