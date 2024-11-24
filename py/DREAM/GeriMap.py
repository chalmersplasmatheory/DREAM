
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


def get(reverse=False, N=256):
    """
    Returns the 'GeriMap' colormap.
    """
    gm = [(0, 0, 0), (.15, .15, .5), (.3, .15, .75),
          (.6, .2, .50), (1, .25, .15), (.9, .5, 0),
          (.9, .75, .1), (.9, .9, .5), (1, 1, 1)]

    if reverse:
        return LinearSegmentedColormap.from_list('GeriMap_r', gm[::-1], N=N)
    else:
        return LinearSegmentedColormap.from_list('GeriMap', gm, N=N)


def register():
    """
    Register the perceptually uniform colormap 'GeriMap' with matplotlib.
    """
    mpl.colormaps.register(cmap=get(False))
    mpl.colormaps.register(cmap=get(True))


