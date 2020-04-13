
import numpy as np
from Equations.EquationException import EquationException


class HotElectronDistribution:
    
    def __init__(self,
        init=None, initr=None, initp=None, initxi=None,
        initppar=None, initpperp=None):
        """
        Constructor.
        """
        self.init = {
            'x': None,
            'r': np.array([]),
            'p': np.array([]),
            'xi': np.array([]),
            'ppar': np.array([]),
            'pperp': np.array([])
        }

        if init is not None:
            self.setInitialValue(init, r=initr, p=initp, xi=initxi, ppar=initppar, pperp=initpperp)


    def setInitialValue(self, init,
        r, p=None, xi=None, ppar=None, pperp=None):
        """
        Set the initial value of this hot electron
        distribution function.

        init:  Array representing the distribution function value on the grid
               (must have size (nr, nxi, np) or (nr, npperp, nppar))
        r:     Radial grid on which the initial distribution is given.
          MOMENTUM GRID SETTINGS:
        p:     Momentum grid.
        xi:    Pitch grid.
          OR
        ppar:  Parallel momentum grid.
        pperp: Perpendicular momentum grid.
        """
        if p is not None and xi is not None:
            self.init['x'] = init
            self.init['r'] = r
            self.init['p'] = p
            self.init['xi'] = xi
            self.init['ppar'] = np.array([])
            self.init['pperp'] = np.array([])
        elif ppar is not None and pperp is not None:
            self.init['x'] = init
            self.init['r'] = r
            self.init['ppar'] = ppar
            self.init['pperp'] = pperp
            self.init['p'] = np.array([])
            self.init['xi'] = np.array([])
        else:
            raise EquationException("f_hot: No momentum grid given for initial value.")

        self.verifyInitialDistribution()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this HotElectronDistribution object.
        """
        data = {'init': {}}

        if self.init['x'] is not None:
            data['init']['x'] = self.init['x']
            data['init']['r'] = self.init['r']

            if self.init['p'].size > 0 and self.init['xi'].size > 0:
                data['init']['p'] = self.init['p']
                data['init']['xi'] = self.init['xi']
            elif self.init['ppar'].size > 0 and self.init['pperp'].size > 0:
                data['init']['ppar'] = self.init['ppar']
                data['init']['pperp'] = self.init['pperp']

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        self.verifyInitialDistribution()


    def verifyInitialDistribution(self):
        """
        Verifies that the initial distribution function has
        been set correctly and consistently.
        """
        if self.init['x'] is None: return

        nr = self.init['r'].size
        p1, p2 = None, None
        p1name, p2name = None, None
        np1, np2 = 0, 0

        if self.init['p'].size > 0 and self.init['xi'].size > 0:
            p1name = 'p'
            p2name = 'xi'
        elif self.init['ppar'].size > 0 and self.init['pperp'].size > 0:
            p1name = 'ppar'
            p2name = 'pperp'
        else:
            raise EquationException("f_hot: No momentum grid given for initial value.")

        p1 = self.init[p1name]
        p2 = self.init[p2name]

        if len(p1.shape) != 1:
            raise EquationException("f_hot: Invalid dimensions of momentum grid '{}'. Must be 1D array.".format(p1name))
        elif len(p2.shape) != 1:
            raise EquationException("f_hot: Invalid dimensions of momentum grid '{}'. Must be 1D array.".format(p2name))

        np1 = p1.size
        np2 = p2.size

        if self.init['x'].shape != (nr, np2, np1):
            raise EquationException("f_hot: Invalid size of initial distribution function: {}. Expected: {}.".format(self.init['x'].shape, (nr, np2, np1)))


