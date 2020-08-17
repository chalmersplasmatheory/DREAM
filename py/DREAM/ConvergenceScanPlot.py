# A class for plotting the result of a
# convergence scan.
############################################

import matplotlib.pyplot as plt
import numpy as np
from . import DREAMIO
from . import GeriMap

from .ConvergenceScan import ConvergenceScan, ConvergenceScanException


class ConvergenceScanPlot:
    
    def __init__(self, scan=None):
        """
        Constructor.

        :param ConvergenceScan scan: ConvergenceScan objec to plot results from.
        """
        self.result = None
        self.outputParameters = None
        
        self.loadResult(scan)


    def loadResult(self, scan):
        """
        Load a ConvergenceScan result to analyse.

        :param scan: Either a ConvergenceScan object or a string. If the former, the result of that scan is loaded. Otherwise, the result is assumed to be stored in the file with the name given by the string.
        """
        if type(scan) is ConvergenceScan:
            self.result = scan.result
            self.outputParameters = scan.getOutputParameters()
        elif type(scan) is str:
            d = DREAMIO.LoadHDF5AsDict(scan)
            self.result = d['result']
            self.outputParameters = d['outputParameters']
        else:
            print("WARNING: Input 'scan' parameter has an unrecognized type. Ignoring...")

        
    def isConverged(self, yesno=True, runIndex=0):
        """
        Checks whether the ConvergenceScan result is converged
        in all output parameters with respect to the scan parameter
        'scanParameter'. If 'scanParameter' is 'None', all scan
        parameters are returned.

        :param bool yesno: If ``True``, returns a single scalar value declaring whether or not all scans are converged.

        :return: If ``yesno`` is ``False``, returns a dict of scan parameters, each of which in turn contains dicts of output parameters, which are either ``True`` or ``False`` depending on if the scan was converged in that output parameter or not.
        :rtype: dict
        """

        conv = {}
        allConverged = True
        for scanParameter in self.result:
            conv[scanParameter] = {}

            for outParameter, op in self.result[scanParameter].items():
                i0 = np.where(op['index'] == runIndex)

                # We prefer to compare to a more well-resolved result
                i1 = None
                try: i1 = np.where(op['index'] == runIndex+1)
                except ValueError: pass

                # If a more well-resolved result is not available, try a
                # less resolved result
                if i1 is None:
                    i1 = np.where(op['index'] == runIndex-1)

                v0 = op['outval'][i0]
                v1 = op['outval'][i1]
                        
                reltol = self.outputParameters[outParameter]['reltol'][0]

                Delta = 1
                if v1 == 0:
                    Delta = np.abs(v0)
                else:
                    Delta = np.abs(v0/v1-1)

                converged = Delta < reltol
                conv[scanParameter][outParameter] = converged

                allConverged = allConverged and converged

        if yesno:
            return allConverged
        else:
            return conv


    def plot(self, plotShape=None, subplot=True, combineOutput=True, normalized=False):
        """
        Plots the result of the convergence scan.

        :param tuple plotShape:    If ``subplot`` is ``True``, sets the shape of the subplot to create (i.e. a tuple with two integer values). The product of the elements of this tuple must be greater than or equal to the number of scan parameters (times the number of output parameters if 'combineOutput' is False).
        :param bool subplot:       If ``True``, plots everything in a single window (but on separate axes).
        :param bool combineOutput: If ``True``, plots the convergence of all output parameters on the same axes. This forces ``normalized = True`` (but only if there are several output parameters).
        :param bool normalized:    If ``True``, plots relative error in the parameter (compared to the most well-resolved run). Otherwise, the output parameters are plotted in absolute units.  (This parameter is automatically forced to True if ``combineOutput`` is ``True``).
        """
        if self.result is None:
            print("WARNING: No result has been loaded.")
            return

        # Verify plot shape
        nscan = len(list(self.result.keys()))
        nout  = len(list(self.outputParameters.keys()))
        m, n  = None, None

        # If there's only one output parameter, 'combineOutput'
        # doesn't make much sense, so to avoid confusing the user
        # with the behaviour of 'normalized' when 'combineOutput = True',
        # we disable 'combineOutput' here.
        if nout == 1:
            combineOutput = False

        colormap = GeriMap.get()

        if plotShape is None:
            if combineOutput:
                m = int(np.ceil(np.sqrt(nscan)))
                n = int(np.ceil(nscan / m))
                plotShape = (m, n)
            else:
                m = nscan
                n = nout
                plotShape = (nout, nscan)
        else:
            x = plotShape[0]*plotShape[1]

            if (combineOutput and x < nscan) or (not combineOutput and x < nscan*nout):
                raise ConvergenceScanException("The specified sub-plot shape is not capable of holding all plots.")

        figsize  = (m*4, n*3)

        # Create figure and axes
        fig, axes = plt.subplots(plotShape[0], plotShape[1], figsize=figsize)
        if plotShape[0] == 1:
            axes = [axes]
        if plotShape[1] == 1:
            axes = [[ax] for ax in axes]

        # Generate plots
        scanParameters = list(self.result.keys())
        outputParameters = list(self.outputParameters.keys())

        for i in range(0, nscan):
            idx = 0
            if combineOutput:
                (I,J) = np.unravel_index(i, plotShape)
                ax = axes[I][J]

                for j in range(0, nout):
                    self.__plotParameter(ax, scanParameters[i], outputParameters[j], normalized=True, ylabel='Relative error', color=colormap(j/(nout+1) * 255))
            else:
                for j in range(0, nout):
                    (I,J) = np.unravel_index(i*nout + j, plotShape)
                    ax = axes[I][J]

                    self.__plotParameter(ax, scanParameters[i], outputParameters[j], normalized=normalized, ylabel=outputParameters[j], color=colormap(0) * 255)

        # Hide remaining axes
        for i in range(nscan*nout, m*n):
            (I, J) = np.unravel_index(i, plotShape)
            fig.delaxes(axes[I][J])

        plt.tight_layout()
        plt.show()


    def __plotParameter(self, ax, scanParameter, outParameter, normalized, ylabel, color):
        """
        Plot a single scan for an individual output parameter.

        ax:            matplotlib axes object to plot on.
        scanParameter: Name of scan parameter to plot.
        outParameter:  Name of output parameter to plot.
        normalized:    Whether or not to plot relative error rather
                       than actual value of parameter.
        """
        i = self.result[scanParameter][outParameter]['index']
        l = [str(x) for x in self.result[scanParameter][outParameter]['scanval']]
        v = self.result[scanParameter][outParameter]['outval']

        if normalized:
            if v[-1] != 0:
                v = (v/v[-1] - 1) * 100
                ylabel += r' (\%)'

        #ax.plot(i, v, linewidth=2, marker='s', markersize=6, color=color)
        ax.plot(i, v, linewidth=2, marker='s', markersize=6)
        ax.set_xlabel(scanParameter)
        ax.set_ylabel(ylabel)
        ax.set_xticks(i)
        ax.set_xticklabels(l)


