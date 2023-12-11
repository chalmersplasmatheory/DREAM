# Container class for linear solver data


import matplotlib.pyplot as plt
import numpy as np
from .Solver import Solver


class SolverNonLinear(Solver):
    

    def __init__(self, solverdata=None, output=None):
        """
        Constructor.
        """
        super().__init__(solverdata, output)

        self.iterations = [int(x) for x in solverdata['iterations'][:]]
        self.backupinverter = [x==1 for x in solverdata['backupinverter'][:]]
        self.nontrivials = solverdata['nontrivials'][:].split(';')[:-1]
        self.unknowns = solverdata['unknowns'][:].split(';')[:-1]

        if 'convergence' in solverdata:
            self.convergence_residual = solverdata['convergence']['residual'][:]


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


    def getBackupRanges(self):
        """
        Return the time step ranges for which the backup solver was used.
        This method returns an array of tuples, where each tuple denotes a
        single range of time steps where the backup solver was used.
        """
        r = np.linspace(1, self.output.grid.t.size, self.output.grid.t.size)[np.where(self.backupinverter)]
        if not r:
            return []

        arr = []
        start = r[0]
        i = 1
        while i < r.size:
            # One step in between _without_ backup solver
            if r[i] > r[i-1]+1:
                arr.append((start, r[i-1]))
                start = r[i]

            i += 1

        arr.append((start, r[-1]))

        return arr

    
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

        # Plot where backup solver is used
        xr = self.getBackupRanges()
        for rg in xr:
            if time:
                ts1 = 0.5*self.output.grid.t[int(rg[0])] + 0.5*self.output.grid.t[int(rg[0])-1]
                ts2 = 0.5*self.output.grid.t[int(rg[1])]

                if rg[1]+1 < self.output.grid.t.size:
                    ts2 += 0.5*self.output.grid.t[int(rg[1])+1]
                else:
                    ts2 += 2*ts2 - 0.5*self.output.grid.t[int(rg[1])-1]
            else:
                ts1, ts2 = rg[0]-0.5, rg[1]+0.5

            ax.fill_between([ts1, ts2], [0, 0], [ymax, ymax], color='r', alpha=0.3)

        if show:
            plt.show(block=False)

        return ax


