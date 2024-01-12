# Class for interacting with solver convergence data

import matplotlib.pyplot as plt
import numpy as np


class SolutionConvergenceInformation:
    

    def __init__(self, x, dx, nontrivials, epsa, epsr):
        """
        Constructor.
        """
        self.x = x
        self.dx = dx
        self.nontrivials = nontrivials

        self.epsa = epsa
        self.epsr = epsr


    def _getSolutionData(self, t, unknown=None):
        """
        Return data for the specified subset.
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
        _x = np.zeros((len(uids), self.x.shape[2]))
        _dx = np.zeros((len(uids), self.dx.shape[2]))
        _epsa = np.zeros((len(uids),))
        _epsr = np.zeros((len(uids),))
        for i in range(len(uids)):
            _x[i,:] = self.x[uids[i],t,:]
            _dx[i,:] = self.dx[uids[i],t,:]
            _epsa[i] = self.epsa[uids[i]]
            _epsr[i] = self.epsr[uids[i]]

        return range(_x.shape[1]), _x, _dx, uids, unknown, _epsa, _epsr


    def plot(self, t, u=None, normalized=True, criterion=True, ax=None, show=None):
        """
        Plot the solution convergence progress.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        iters, x, dx, uids, unknowns, epsa, epsr = self._getSolutionData(t=t, unknown=u)

        if normalized:
            dxx = dx / x
        else:
            dxx = dx

        for i in range(len(uids)):
            p = ax.semilogy(iters, dxx[i,:], '.-', label=unknowns[i])

            if criterion:
                if normalized:
                    epsax = epsa[i] / x[i,:]
                    epsrx = epsr[i] * np.ones(len(iters))
                else:
                    epsax = np.ones(len(iters)) * epsa[i]
                    epsrx = epsr[i] * x[i,:]

                ax.semilogy(iters, epsax + epsrx, '--', color=p[0].get_color())

        ax.set_xlabel('Iterations')
        if normalized:
            ax.set_ylabel('$|dx| / |x|$')
        else:
            ax.set_ylabel('$|dx|$')

        ax.legend(fancybox=False, edgecolor='k')

        if show:
            plt.show()

        return ax


    def plotx(self, t, u=None, ax=None, show=None):
        """
        Plot the solution convergence progress.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        iters, x, _, uids, unknowns, epsa, epsr = self._getSolutionData(t=t, unknown=u)

        if normalized:
            xx = x - x[-1]
        else:
            xx = x

        for i in range(len(uids)):
            ax.semilogy(iters, xx[i,:], '.-', label=unknowns[i])

        ax.set_xlabel('Iterations')
        if normalized:
            ax.set_ylabel(r'$|x-x^\star|$')
        else:
            ax.set_ylabel('$|x|$')

        ax.legend(fancybox=False, edgecolor='k')

        if show:
            plt.show()

        return ax


