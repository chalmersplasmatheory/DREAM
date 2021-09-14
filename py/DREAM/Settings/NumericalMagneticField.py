# Base class for numerical magnetic fields


class NumericalMagneticField:
    

    def __init__(self, a, R0):
        """
        Constructor.
        """
        self.a = a
        self.R0 = R0

    
    def visualize(self, ax=None, show=None):
        """
        Visualize this numerical magnetic field.
        """
        raise Exception("This numerical magnetic field does not implement 'visualize()'.")


