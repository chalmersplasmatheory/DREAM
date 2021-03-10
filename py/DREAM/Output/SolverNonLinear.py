# Container class for linear solver data


from .Solver import Solver


class SolverNonLinear(Solver):
    

    def __init__(self, solverdata=None, output=None):
        """
        Constructor.
        """
        super().__init__(solverdata, output)

        self.iterations = [int(x) for x in solverdata['iterations'][:]]
        self.backupinverter = [x==1 for x in solverdata['backupinverter'][:]]


    def __str__(self):
        """
        Convert this object into a string.
        """
        s = "Non-linear solver statistics\n\n"

        s += "Max. iterations: {}\n".format(max(self.iterations))
        s += "Avg. iterations: {}\n".format(sum(self.iterations)/len(self.iterations))
        s += "Min. iterations: {}\n\n".format(min(self.iterations))
        
        bi = sum(self.backupinverter)
        if bi == 0:
            s += "Backup inverter not used\n"
        elif bi == 1:
            s += "Backup inverter used: 1 time\n"
        else:
            s += "Backup inverter used: {} times\n".format(bi)

        return s

    
    def plot(self, ax=None, show=None, **kwargs):
        """
        Visualize the solver statistics.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()
            
            if show is None:
                show = True

        ax.plot(self.output.grid.t[1:], self.iterations, linewidth=2, **kwargs)
        ax.set_xlabel(r'Simulation time (s)')
        ax.set_ylabel(r'Number of iterations')

        if show:
            plt.show(block=False)

        return ax


