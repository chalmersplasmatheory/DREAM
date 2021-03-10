# Container class for linear solver data


from .Solver import Solver


class SolverLinear(Solver):
    

    def __init__(self, solverdata=None, output=None):
        """
        Constructor.
        """
        super().__init__(solverdata, output)


    def __str__(self):
        """
        Convert this object into a string.
        """
        return "<Empty linear solver data>"


