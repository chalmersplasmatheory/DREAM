# Container class for linear solver data


import matplotlib.pyplot as plt
from .Solver import Solver


class SolverLinear(Solver):
    

    def __init__(self, solverdata=None, output=None):
        """
        Constructor.
        """
        super().__init__(solverdata, output)

        if 'iterations' in solverdata:
            self.iterations = [int(x) for x in solverdata['iterations'][:]]
        else:
            self.iterations = [1 for _ in range(output.grid.t.size)]


    def __str__(self):
        """
        Convert this object into a string.
        """
        return "<Empty linear solver data>"


    def plot(self, time=True, ax=None, show=None, **kwargs):
        """
        Visualize the solver statistics.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()
            
            if show is None:
                show = True

        if time:
            t = self.output.grid.t[1:]
            xlbl = r'Simulation time (s)'
        else:
            t = np.linspace(1, self.output.grid.t.size-1, self.output.grid.t.size-1)
            xlbl = r'Time step'

        ax.plot(t, self.iterations, linewidth=2, **kwargs)

        ax.set_xlabel(xlbl)
        ax.set_ylabel(r'Number of iterations')

        ymax = max(self.iterations)*1.2
        ax.set_ylim([0, ymax])

        if show:
            plt.show(block=False)

        return ax



