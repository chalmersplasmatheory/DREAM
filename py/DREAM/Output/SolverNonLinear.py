# Container class for linear solver data


import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
from .Solver import Solver
from .SolutionConvergenceInformation import SolutionConvergenceInformation


class SolverNonLinear(Solver):
    

    def __init__(self, solverdata=None, output=None):
        """
        Constructor.
        """
        super().__init__(solverdata, output)

        if 'solvertime' in solverdata:
            self.solvertime = solverdata['solvertime'][:]
        if 'iterations' in solverdata:
            self.iterations = solverdata['iterations'][:].astype(int)
        if 'backupinverter' in solverdata:
            self.backupinverter = solverdata['backupinverter'][:].astype(bool)
        if 'nontrivials' in solverdata:
            self.nontrivials = solverdata['nontrivials'][:].split(';')[:-1]
        if 'unknowns' in solverdata:
            self.unknowns = solverdata['unknowns'][:].split(';')[:-1]

        if 'convergence' in solverdata:
            self.convergence_residual = solverdata['convergence']['residual'][:]
            self.convergence_residualmaxerr = solverdata['convergence']['residualmaxerror'][:]

            if 'x' in solverdata['convergence']:
                self.solution = SolutionConvergenceInformation(
                    solverdata['convergence']['x'][:],
                    solverdata['convergence']['dx'][:],
                    self.nontrivials,
                    epsa=solverdata['convergence']['epsa'][:],
                    epsr=solverdata['convergence']['epsr'][:]
                )


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
        r = np.linspace(1, self.solvertime.size, self.solvertime.size)[np.where(self.backupinverter)]
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
            t = self.solvertime[:]
            xlbl = r'Simulation time (s)'
        else:
            t = np.linspace(1, self.solvertime.size, self.solvertime.size)
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
                ts1 = 0.5*self.solvertime[int(rg[0])] + 0.5*self.solvertime[int(rg[0])-1]
                ts2 = 0.5*self.solvertime[int(rg[1])]

                if rg[1]+1 < self.solvertime.size:
                    ts2 += 0.5*self.solvertime[int(rg[1])+1]
                else:
                    ts2 += 2*ts2 - 0.5*self.solvertime[int(rg[1])-1]
            else:
                ts1, ts2 = rg[0]-0.5, rg[1]+0.5

            ax.fill_between([ts1, ts2], [0, 0], [ymax, ymax], color='r', alpha=0.3)

        if show:
            plt.show(block=False)

        return ax


    def _getResidualData(self, data, unknown=None, t=None, showtimesteps=False):
        """
        Returns the requested subset of residual data.
        """
        if unknown is None:
            unknown = self.nontrivials
            uids = range(len(self.nontrivials))
        else:
            if type(unknown) != list and type(unknown) != tuple:
                unknown = [unknown]

            # Translate unknown names to indices
            uids = []
            for u in unknown:
                uids.append(self.nontrivials.index(u))

        # Select data to plot
        if t is not None:
            _d = np.zeros((len(uids), 1))
            if showtimesteps:
                tarr = np.array([t])
            else:
                tarr = np.array([self.solvertime[t-1]])

            for i in range(len(uids)):
                _d[i,0] = data[uids[i],t]
        else:
            _d = np.zeros((len(uids), data.shape[1]))
            if showtimesteps:
                tarr = np.array(range(self.solvertime.size))+1
            else:
                tarr = self.solvertime[:]
            for i in range(len(uids)):
                _d[i,:] = data[uids[i],:]

        # Plot
        if _d.shape[0] == 1:
            idxs = [-1, 1]
            _d = np.matlib.repmat(_d, 2, 1)
        else:
            idxs = list(range(len(uids)))

        return tarr, idxs, _d, uids, unknown

        
    def plotResidualConvergence(self, unknown=None, t=None, ax=None, show=None, showtimesteps=False):
        """
        Plot the residual convergence status for one or more unknowns.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()
            
            if show is None:
                show = True

        tarr, idxs, data, uids, unknowns = self._getResidualData(self.convergence_residual, unknown=unknown, t=t, showtimesteps=showtimesteps)
        ax.pcolormesh(tarr, idxs, data, cmap='RdYlGn', shading='nearest', vmin=0, vmax=1)
        ax.set_yticks(list(range(len(uids))), labels=unknowns)

        if show:
            plt.show(block=False)

        return ax


    def plotResidualMaxError(
        self, unknown=None, t=None, ax=None, show=None,
        showtimesteps=False, logx=False
    ):
        """
        Plot the maximum error in the residual.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        tarr, _, data, uids, unknowns = self._getResidualData(self.convergence_residualmaxerr, unknown=unknown, t=t, showtimesteps=showtimesteps)

        ax.plot([tarr[0], tarr[-1]], [1, 1], 'k--')
        for i in range(len(unknowns)):
            ax.plot(tarr, data[i,:], label=unknowns[i])

        ax.legend(fancybox=False, edgecolor='k')
        if logx:
            ax.set_xscale('log')
        ax.set_yscale('log')

        if show:
            plt.show()

        return ax


