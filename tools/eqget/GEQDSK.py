# 
# Class for loading and working with a magnetic equilibrium stored in a
# GEQDSK file.
#
# Written by: Mathias Hoppe, 2022
#

import h5py
from matplotlib._contour import QuadContourGenerator
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.interpolate import CubicSpline, InterpolatedUnivariateSpline, RectBivariateSpline
from EqBase import EqBase


class GEQDSK(EqBase):
    

    def __init__(self, filename, cocos=1):
        """
        Constructor.
        """
        self.load(filename, cocos=cocos)


    def load(self, filename, cocos=1):
        """
        Load data from the named GEQDSK file to this object.
        """
        data = self.load_geqdsk(filename, cocos=1)
        self.process_data(data)


    def load_geqdsk(self, filename, cocos=1):
        """
        Load the named GEQDSK file.
        """
        with open(filename) as fh:
            header = fh.readline()
            words = header.split()
            if len(words) < 3:
                raise ValueError("Expected at least 3 numbers on first line")

            nx, ny = int(words[-2]), int(words[-1])
            
            data = {"nx": nx, "ny": ny, "cocos": cocos}
            fields = ["rboxlen", "zboxlen", "rcentr", "rleft", "zmid", "raxis",
                      "zaxis", "psiaxis", "psiedge", "bcentr", "cpasma", "psiaxis",
                      None, "raxis", None, "zaxis", None, "psiedge", None, None]

            values = self._next_value(fh)
            
            for f in fields:
                val = next(values)
                if f:
                    data[f] = val

            def _read_1d(n):
                """
                Read a 1D array of length n from the GEQDSK file.
                """
                val = np.zeros(n)
                for i in range(n):
                    val[i] = next(values)

                return val


            def _read_2d(n, m):
                """
                Read a 2D (n,m) array in Fortran order
                """
                val = np.zeros((n, m))
                for j in range(m):
                    for i in range(n):
                        val[i, j] = next(values)

                return val


            data["fpol"] = _read_1d(nx)
            data["p"] = _read_1d(nx)
            data["ffprime"] = _read_1d(nx)
            data["pprime"] = _read_1d(nx)

            data["psi"] = _read_2d(nx, ny)

            data["q"] = _read_1d(nx)

            # Ensure that psi is divided by 2pi
            if cocos > 10:
                for var in ["psi", "psiaxis", "psiedge"]:
                    data[var] /= 2 * pi

            nbdry = next(values)
            nlim = next(values)

            if nbdry > 0:
                data["rbdry"] = np.zeros(nbdry)
                data["zbdry"] = np.zeros(nbdry)
                for i in range(nbdry):
                    data["rbdry"][i] = next(values)
                    data["zbdry"][i] = next(values)

            if nlim > 0:
                data["rlim"] = np.zeros(nlim)
                data["zlim"] = np.zeros(nlim)
                for i in range(nlim):
                    data["rlim"][i] = next(values)
                    data["zlim"][i] = next(values)

            return data


