# A class for generating matplotlib plots from DREAMOutput objects.

import matplotlib.pyplot as plt
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
        
        self.config = {}


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

        self.config.append({
            'name': name,
            'x': x,
            'y': y
        })

        if title is not None: config['title'] = title
        if xlabel is not None: config['xlabel'] = xlabel
        if xmin is not None: config['xmin'] = xmin
        if xmax is not None: config['xmax'] = xmax
        if ymin is not None: config['ymin'] = ymin
        if ymax is not None: config['ymax'] = ymax
        if log is not None: config['log'] = log


    def removeQuantity(self, name):
        """
        Remove a quantity from the plot configuration.
        """
        # Locate item to delete
        found = False
        for i in range(len(self.config)):
            if self.config[i]['name'] == name:
                del self.config[i]
                break

        if not found:
            raise Exception(f"The quantity with name '{name}' has not been configured.")


    def render(self, output):
        """
        Render the plot with data from the given data object.
        
        :param output: Data to render (can be a dict or a DREAMOutput object).
        """
        isUpdate = True
        if self.axs is None or self.axs.shape != self.shape:
            self.axs = fig.subplots(nrows=self.shape[0], ncols=self.shape[1]).flatten()
            isUpdate = False

        for i in range(len(self.config)):
            self._makeplot(output, self.axs[i], self.config[i], isSimpleUpdate=isUpdate)


    def setShape(self, nrows, ncols):
        """
        Set the number of rows and columns to use for the figure.
        """
        if not isinstance(nrows, numbers.Number):
            raise Exception("The number of rows 'nrows' must be a number.")
        if not isinstance(ncols, numbers.Number):
            raise Exception("The number of columns 'ncols' must be a number.")

        self.shape = (nrows, ncols)
    

    def _makeplot(self, output, ax, config, isSimpleUpdate):
        """
        Paint a single plot on the given axis, using the given
        plot configuration.

        :param output: Output object to fetch data from.
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
            x = xf(output)
        elif type(xf) == str:
            # TODO
            pass
        else:
            raise Exception(f"{name}: Unrecognized type of 'x' specification.")

        if callable(yf):
            y = yf(output)
        elif type(yf) == str:
            # TODO
            pass
        else:
            raise Exception(f"{name}: Unrecognized type of 'y' specification.")

        if name in self.lines and isSimpleUpdate:
            h = self.lines[name]
            h.set_data(x, y)
        else:
            if log:
                h, _ = ax.semilogy(x, y)
            else:
                h, _ = ax.plot(x, y)

            if 'title' in config:
                ax.set_title(config['title'])
            if xlabel is not None:
                ax.set_xlabel(xlabel)

            if (xmin is not None) and (xmax is not None):
                ax.set_xlim((xmin, xmax))
            if (ymin is not None) and (ymax is not None):
                ax.set_ylim((ymin, ymax))

            self.lines[name] = h


