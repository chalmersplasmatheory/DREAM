#
# RadialGrid settings object.
########################################

import numpy as np
from DREAM.DREAMException import DREAMException


TYPE_CYLINDRICAL = 1
TYPE_ANALYTIC_TOROIDAL = 2


class RadialGrid:
    

    def __init__(self, ttype=1):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        # Cylindrical settings
        self.a  = 0.0
        self.B0 = 0.0
        self.nr = int(0)
        self.r0 = 0.0

        # Ripple parameters
        self.ripple_ncoils = 0
        self.ripple_deltacoils = 0.0
        self.ripple_m = None
        self.ripple_n = None
        self.ripple_dB_B = None
        self.ripple_r = None
        self.ripple_t = None


    #######################
    # SETTERS
    #######################
    def setB0(self, B0):
        if B0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'B0'.")
        
        self.B0 = float(B0)


    def setInnerRadius(self, r0):
        if r0 < 0:
            raise DREAMException("RadialGrid: Invalid value assigned to innermost radius 'r0': {}".format(r0))

        self.r0 = r0


    def setMinorRadius(self, a):
        if a <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(a))

        self.a = float(a)


    def setNr(self, nr):
        if nr <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'nr': {}".format(nr))

        self.nr = int(nr)


    def setRipple(self, m, n, dB_B, ncoils=0, deltacoils=0, r=[0], t=[0]):
        """
        Enable the ripple pitch scattering term.

        :param list m:           Poloidal mode numbers of magnetic perturbation(s).
        :param list n:           Toroidal mode numbers of magnetic perturbation(s).
        :param dB_B:             Magnetic perturbations (shape: nModes x nt x nr).
        :param int ncoils:       Number of toroidal field coils.
        :param float deltacoils: Distance between toroidal field coils.
        :param r:                Radial grid on which the magnetic perturbations are given.
        :param t:                Time grid on which the magnetic perturbations are given.
        """
        if type(m) == list: m = np.array(m)
        elif np.isscalar(m): m = np.array([float(m)])

        if type(n) == list: n = np.array(n)
        elif np.isscalar(n): n = np.array([float(n)])

        if type(r) == list: r = np.array(r)
        elif np.isscalar(r): r = np.array([float(r)])

        if type(t) == list: t = np.array(t)
        elif np.isscalar(t): t = np.array([float(t)])

        if type(dB_B) == list:
            dB_B = np.array(dB_B)
        if type(dB_B) == float or dB_B.ndim == 1:
            dB_B = np.ones((m.size, t.size, r.size)) * dB_B

        if m.size != n.size:
            raise EquationException("{}: m and n must have the same number of elements.".format(self.name))
        elif dB_B.ndim == 1 and dB_B.size == m.size:
            dB_B = dB_B*np.ones((m.size, t.size, r.size))
        elif dB_B.ndim != 3 or dB_B.shape != (m.size, t.size, r.size):
            raise EquationException("{}: Invalid dimensions of parameter 'dB_B'. Expected {}, but array has {}.".format(self.name, (m.size, t.size, r.size), dB_B.shape))

        self.ripple_ncoils = int(ncoils)
        self.ripple_deltacoils = float(deltacoils)
        self.ripple_m = m
        self.ripple_n = n
        self.ripple_dB_B = dB_B
        self.ripple_r = r
        self.ripple_t = t


    def setType(self, ttype):
        if ttype == TYPE_CYLINDRICAL:
            self.type = ttype
        elif ttype == TYPE_ANALYTIC_TOROIDAL:
            #self.type = ttype
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(ttype))

    
    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.type = data['type']

        if self.type == TYPE_CYLINDRICAL:
            self.a = data['a']
            self.B0 = data['B0']
            self.nr = data['nr']
            self.r0 = data['r0']
        elif self.type == TYPE_ANALYTICAL_TOROIDAL:
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        if 'ripple' in data:
            self.ripple_ncoils = int(scal(data['ripple']['ncoils']))
            self.ripple_deltacoils = float(scal(data['ripple']['deltacoils']))
            self.ripple_m = data['ripple']['m']
            self.ripple_n = data['ripple']['n']
            self.ripple_dB_B = data['ripple']['x']
            self.ripple_r = data['ripple']['r']
            self.ripple_t = data['ripple']['t']


    def todict(self, verify=True):
        """
        Returns the settings in this object as a Python dictionary.
        """
        if verify:
            self.verifySettings()

        data = {
            'type': self.type
        }

        if self.type == TYPE_CYLINDRICAL:
            data['a'] = self.a
            data['B0'] = self.B0
            data['nr'] = self.nr
            data['r0'] = self.r0
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        if self.ripple_ncoils > 0 or self.ripple_deltacoils > 0:
            data['ripple'] = {
                'ncoils': self.ripple_ncoils,
                'deltacoils': self.ripple_deltacoils,
                'm': self.ripple_m,
                'n': self.ripple_n,
                'x': self.ripple_dB_B,
                'r': self.ripple_r,
                't': self.ripple_t
            }

        return data
        
            
    def verifySettings(self):
        """
        Verfiy that the RadialGrid settings are consistent.
        """
        if self.type == TYPE_CYLINDRICAL:
            if self.a is None or self.a <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(self.a))
            elif self.B0 is None or self.B0 <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to 'B0': {}".format(self.B0))
            elif self.r0 is None or self.r0 < 0:
                raise DREAMException("RadialGrid: Invalid value assigned to innermost simulated radius 'r0': {}".format(self.r0))

            if self.nr <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned 'nr': {}. Must be > 0.".format(self.nr))

            if self.r0 >= self.a:
                raise DREAMException("RadialGrid: 'r0' must be strictly less than 'a'.")
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            raise DREAMException("RadialGrid: The analytical toroidal grid has not been implemented yet.")
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        # Ripple settings
        if self.ripple_ncoils > 0 or self.ripple_deltacoils > 0:
            if type(self.ripple_m) != np.ndarray or self.ripple_m.ndim != 1:
                raise EquationException("{}: Invalid type or shape of 'ripple_m'.".format(self.name))
            elif type(self.ripple_n) != np.ndarray or self.ripple_n.ndim != 1:
                raise EquationException("{}: Invalid type or shape of 'ripple_n'.".format(self.name))
            elif self.ripple_m.size != self.ripple_n.size:
                raise EquationException("{}: 'ripple_m' and 'ripple_n' must have the same number of elements.".format(self.name))
            elif type(self.ripple_r) != np.ndarray or self.ripple_r.ndim != 1:
                raise EquationException("{}: Invalid type or shape of 'ripple_r'.".format(self.name))
            elif type(self.ripple_t) != np.ndarray or self.ripple_t.ndim != 1:
                raise EquationException("{}: Invalid type or shape of 'ripple_t'.".format(self.name))
            elif type(self.ripple_dB_B) != np.ndarray or self.ripple_dB_B.shape != (self.ripple_m.size, self.ripple_r.size, self.ripple_t.size):
                raise EquationException("{}: Invalid type or shape of 'ripple_dB_B'.".format(self.ripple_dB_B))
        

