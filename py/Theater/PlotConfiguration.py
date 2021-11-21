# A class for generating matplotlib plots from DataProvider objects.

import matplotlib.pyplot as plt
import numpy as np
import numbers


class PlotConfiguration:
    

    def __init__(self, fig, nrows=1, ncols=1):
        """
        Constructor.
        """
        self.fig = fig
        self.axs = None
        self.lines = {}
        self.shape = (nrows, ncols)

        self.cache = {}
        
        self.config = []


    def addQuantity(self, name, x, y, title=None, xlabel=None, xmin=None, xmax=None,
                    ymin=None, ymax=None, log=False):
        """
        Add a quantity to plot.

        :param name: Name of plot.
        :param x: One of (i) string naming the x variable to use, (ii) string that can be evaluated by the Python interpreter, (iii) Python function
        :param y: One of (i) string that can be evaluated by the Python interpreter, (ii) Python function
        :param title: Figure title to use.
        :param xlabel: Label to use for x axis.
        :param xmin: Lower limit of x axis.
        :param xmax: Upper limit of x axis.
        :param ymin: Lower limit of y axis.
        :param ymax: Upper limit of y axis.
        :param log: If ``True``, use a logarithmic y scale.
        """
        if not callable(x) and type(x) != str:
            raise Exception("The variable 'x' must be either a string or a callable.")
        elif not callable(y) and type(y) != str:
            raise Exception("The variable 'y' must be either a string or a callable.")

        config = {
            'name': name,
            'x': x,
            'y': y
        }

        if title is not None: config['title'] = title
        if xlabel is not None: config['xlabel'] = xlabel
        if xmin is not None: config['xmin'] = xmin
        if xmax is not None: config['xmax'] = xmax
        if ymin is not None: config['ymin'] = ymin
        if ymax is not None: config['ymax'] = ymax
        if log is not None: config['log'] = log

        self.config.append(config)

        # Need to update shape of subplots?
        self.reshape()


    def loadData(self, data):
        """
        Loads the most recently calculated data from the given
        DataProvider.
        """
        self.cache = {
            't': data.getTime(),
            'r': data.getRadius()
        }

        return self.cache


    def monitorsQuantity(self, name):
        """
        Checks if the named quantity has been added to this configuration
        and returns ``True`` if this is the case. Otherwise, returns ``False``.
        """
        for i in range(len(self.config)):
            if self.config[i]['name'] == name:
                return True

        return False


    def removeByAxes(self, ax):
        """
        Remove a quantity by specifying the axes on which it is plotted.
        """
        found = False
        for i in range(len(self.config)):
            if self.axs[i] == ax:
                del self.config[i]
                found = True
                break

        if not found:
            raise Exception(f"The selected quantity could not be found in the plot.")

        self.reshape()


    def removeQuantity(self, name):
        """
        Remove a quantity from the plot configuration.
        """
        # Locate item to delete
        found = False
        for i in range(len(self.config)):
            if self.config[i]['name'] == name:
                del self.config[i]
                found = True
                break

        if not found:
            raise Exception(f"The quantity with name '{name}' has not been configured.")

        self.reshape()


    def render(self, data, clearAxes=False):
        """
        Render the plot with data from the given data object.
        
        :param data: Data to render (should be a DataProvider object).
        :param clearAxes: If ``True``, clears and recreates the figure axes.
        """
        # Pre-load calculated data
        self.loadData(data)

        # Create subplot axes (if needed)
        isUpdate = True
        if self.axs is None or self.axs.shape != self.shape:
            self.fig.clear()
            self.axs = self.fig.subplots(nrows=self.shape[0], ncols=self.shape[1])
            if type(self.axs) != np.array:
                self.axs = np.array(self.axs)
            self.axs = self.axs.flatten()

            isUpdate = False

        # Draw on the axes
        for i in range(len(self.config)):
            self._makeplot(self.axs[i], self.config[i], isSimpleUpdate=isUpdate)

        # Clear empty axes
        if not isUpdate:
            for i in range(len(self.config), self.shape[0]*self.shape[1]):
                self.axs[i].axis('off')

        self.fig.canvas.draw()


    def reshape(self):
        """
        If necessary, reshapes the subplot configuration so that it
        has as few empty plots as possible.
        """
        #if self.shape[0]*self.shape[1] < len(self.config):
        l = len(self.config)
        m = int(np.floor(np.sqrt(l)))
        n = m

        if m*n >= l:
            self.shape = (m, n)
        elif m*(n+1) >= l:
            self.shape = (m, n+1)
        else:
            self.shape = (m+1, n+1)


    def setShape(self, nrows, ncols):
        """
        Set the number of rows and columns to use for the figure.
        """
        if not isinstance(nrows, numbers.Number):
            raise Exception("The number of rows 'nrows' must be a number.")
        if not isinstance(ncols, numbers.Number):
            raise Exception("The number of columns 'ncols' must be a number.")

        self.shape = (nrows, ncols)
    

    def _evalExpression(self, expr):
        """
        Evaluate the given expression.
        """
        local = {
            't': self.cache['t'],
            'r': self.cache['r']
        }

        return eval(expr, globals(), local)


    def _makeplot(self, ax, config, isSimpleUpdate):
        """
        Paint a single plot on the given axis, using the given
        plot configuration.

        :param ax: Axis to plot on.
        :param config: Instructions for how to paint plot.
        :param bool isSimpleUpdate: If ``True``, indicates that this is a minor update of the plot so that 'line.set_data()' is admissible to use. Otherwise the whole plot is completely repainted.
        """
        name   = config['name']
        xf     = config['x']
        yf     = config['y']
        xlabel = None if 'xlabel' not in config else config['xlabel']
        xmin   = None if 'xmin' not in config else config['xmin']
        xmax   = None if 'xmax' not in config else config['xmax']
        ymin   = None if 'ymin' not in config else config['ymin']
        ymax   = None if 'ymax' not in config else config['ymax']
        log    = False if 'log' not in config else config['log']

        if callable(xf):
            x = xf(self.cache[name])
        elif type(xf) == str:
            x = self._evalExpression(xf)
        else:
            raise Exception(f"{name}: Unrecognized type of 'x' specification.")

        if callable(yf):
            y = yf(self.cache[name])
        elif type(yf) == str:
            y = self._evalExpression(xf)
        else:
            raise Exception(f"{name}: Unrecognized type of 'y' specification.")

        if name in self.lines and isSimpleUpdate:
            h = self.lines[name]
            h.set_data(x, y)
        else:
            if log:
                h = ax.semilogy(x, y)
            else:
                h = ax.plot(x, y)

            if 'title' in config:
                ax.set_title(config['title'])
            else:
                ax.set_title(name)
            if xlabel is not None:
                ax.set_xlabel(xlabel)

            if (xmin is not None) and (xmax is not None):
                ax.set_xlim((xmin, xmax))
            if (ymin is not None) and (ymax is not None):
                ax.set_ylim((ymin, ymax))

            self.lines[name] = h


