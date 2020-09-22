# General transport settings which can be added to any UnknownQuantity


import numpy as np
from .. DREAMException import DREAMException


TRANSPORT_NONE = 1
TRANSPORT_PRESCRIBED = 2

BC_CONSERVATIVE = 1
BC_F_0 = 2


class TransportSettings:
    

    def __init__(self, kinetic=False):
        """
        Constructor.

        :param bool kinetic: If ``True``, the coefficient will be assumed kinetic (4D). Otherwise fluid (2D).
        """
        self.kinetic = kinetic
        self.type    = TRANSPORT_NONE

        self.ar        = None
        self.ar_t      = None
        self.ar_r      = None
        self.ar_p      = None
        self.ar_xi     = None
        self.ar_ppar   = None
        self.ar_pperp  = None

        self.drr       = None
        self.drr_t     = None
        self.drr_r     = None
        self.drr_p     = None
        self.drr_xi    = None
        self.drr_ppar  = None
        self.drr_pperp = None

        self.boundarycondition = BC_CONSERVATIVE


    def isKinetic(self): return self.kinetic
    

    def prescribeAdvection(self, ar, t=None, r=None, p=None, xi=None, ppar=None, pperp=None):
        """
        Set the advection coefficient to use.
        """
        self._prescribeCoefficient('ar', coeff=ar, t=t, r=r, p=p, xi=xi, ppar=ppar, pperp=pperp)


    def prescribeDiffusion(self, drr, t=None, r=None, p=None, xi=None, ppar=None, pperp=None):
        """
        Set the diffusion coefficient to use.
        """
        self._prescribeCoefficient('drr', coeff=drr, t=None, r=None, p=None, xi=None, ppar=None, pperp=None)


    def _prescribeCoefficient(self, name, coeff, t=None, r=None, p=None, xi=None, ppar=None, pperp=None):
        """
        General method for prescribing an advection or diffusion coefficient.
        """
        self.type = TRANSPORT_PRESCRIBED

        if np.isscalar(coeff):
            r = np.array([0])
            t = np.array([0])
            p = np.array([0])
            xi = np.array([0])

            if self.kinetic:
                coeff = coeff * np.ones((1,)*4)
            else:
                coeff = coeff * np.ones((1,)*2)

        r = np.asarray(r)
        t = np.asarray(t)

        if self.kinetic == False and len(coeff.shape) == 2:
            setattr(self, name, coeff)
            setattr(self, name+'_r', r)
            setattr(self, name+'_t', t)
        elif self.kinetic == True and len(coeff.shape) == 4:
            # Verify that the momentum grid is given
            if p is not None and xi is not None:
                ppar, pperp = None, None
            elif ppar is not None and pperp is not None:
                p, xi = None, None
            else:
                raise TransportException("No momentum grid provided for the 4D transport coefficient.")

            setattr(self, name, coeff)
            setattr(self, name+'_r', r)
            setattr(self, name+'_t', t)

            if p is not None:
                setattr(self, name+'_p', r)
                setattr(self, name+'_xi', r)
            else:
                setattr(self, name+'_ppar', ppar)
                setattr(self, name+'_pperp', pperp)
        else:
            raise TransportException("Invalid dimensions of prescribed coefficient: {}. Expected {} dimensions.".format(coeff.shape, 4 if self.kinetic else 2))


    def setBoundaryCondition(self, bc):
        """
        Set the boundary condition to use for the transport.
        """
        self.boundarycondition = bc


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        self.ar = None
        self.ar_r = None
        self.ar_t = None
        self.ar_p = None
        self.ar_xi = None
        self.ar_ppar = None
        self.ar_pperp = None

        self.drr = None
        self.drr_r = None
        self.drr_t = None
        self.drr_p = None
        self.drr_xi = None
        self.drr_ppar = None
        self.drr_pperp = None

        if 'type' in data:
            self.type = data['type']

        if 'boundarycondition' in data:
            self.boundarycondition = data['boundarycondition']

        if 'ar' in data:
            self.ar = data['ar']['x']
            self.r  = data['ar']['r']
            self.t  = data['ar']['t']

            if self.kinetic:
                if 'p' in data['ar']: self.ar_p = data['ar']['p']
                if 'xi' in data['ar']: self.ar_xi = data['ar']['xi']
                if 'ppar' in data['ar']: self.ar_ppar = data['ar']['ppar']
                if 'pperp' in data['ar']: self.ar_pperp = data['ar']['pperp']

        if 'drr' in data:
            self.drr = data['drr']['x']
            self.drr_r  = data['drr']['r']
            self.drr_t  = data['drr']['t']

            if self.kinetic:
                if 'p' in data['drr']: self.drr_p = data['drr']['p']
                if 'xi' in data['drr']: self.drr_xi = data['drr']['xi']
                if 'ppar' in data['drr']: self.drr_ppar = data['drr']['ppar']
                if 'pperp' in data['drr']: self.drr_pperp = data['drr']['pperp']


    def todict(self):
        """
        Returns these settings as a dictionary.
        """
        data = {
            'type': self.type,
            'boundarycondition': self.boundarycondition
        }

        # Advection?
        if self.type == TRANSPORT_PRESCRIBED and self.ar is not None:
            data['ar'] = {
                'x': self.ar,
                'r': self.ar_r,
                't': self.ar_t
            }

            if self.kinetic:
                if self.ar_p is not None:
                    data['ar']['p'] = self.ar_p
                    data['ar']['xi'] = self.ar_xi
                else:
                    data['ar']['ppar'] = self.ar_ppar
                    data['ar']['pperp'] = self.ar_pperp

        # Diffusion?
        if self.type == TRANSPORT_PRESCRIBED and self.drr is not None:
            data['drr'] = {
                'x': self.drr,
                'r': self.drr_r,
                't': self.drr_t
            }

            if self.kinetic:
                if self.drr_p is not None:
                    data['drr']['p'] = self.drr_p
                    data['drr']['xi'] = self.drr_xi
                else:
                    data['drr']['ppar'] = self.drr_ppar
                    data['drr']['pperp'] = self.drr_pperp
        
        return data


    def verifySettings(self):
        """
        Verify that the settings are consistent.
        """
        if self.type == TRANSPORT_NONE:
            pass
        elif self.type == TRANSPORT_PRESCRIBED:
            self.verifySettingsPrescribedCoefficient('ar')
            self.verifySettingsPrescribedCoefficient('drr')

            bcs = [BC_CONSERVATIVE, BC_F_0]
            if self.boundarycondition not in bcs:
                raise TransportException("{}: Invalid boundary condition specified for transport: {}".format(coeff, self.boundarycondition))
        else:
            raise TransportException("Unrecognized transport type: {}".format(self.type))


    def verifySettingsCoefficient(self, coeff):
        """
        Verify consistency of the named prescribed transport coefficient.
        """
        g = lambda v : self.__dict__[coeff+v]
        c = g('')

        if self.kinetic:
            if len(c.shape) != 4:
                raise TransportException("{}: Invalid dimensions of transport coefficient: {}".format(coeff, c.shape))
            elif not np.isarray(g('_t')) or len(g('_t').shape) != 1 or g('_t').size != c.shape[0]:
                raise TransportException("{}: Invalid dimensions of time vector. Expected {} elements.".format(coeff, c.shape[0]))
            elif not np.isarray(g('_r')) or len(g('_r').shape) != 1 or g('_r').size != c.shape[1]:
                raise TransportException("{}: Invalid dimensions of radius vector. Expected {} elements.".format(coeff, c.shape[1]))

            if g('_p') is not None or g('_xi') is not None:
                if not np.isarray(g('_xi')) or len(g('_xi').shape) != 1 or g('_xi').size != c.shape[2]:
                    raise TransportException("{}: Invalid dimensions of xi vector. Expected {} elements.".format(coeff, c.shape[2]))
                elif not np.isarray(g('_p')) or len(g('_p').shape) != 1 or g('_p').size != c.shape[3]:
                    raise TransportException("{}: Invalid dimensions of p vector. Expected {} elements.".format(coeff, c.shape[3]))
            elif g('_ppar') is not None or g('_pperp') is not None:
                if not np.isarray(g('_pperp')) or len(g('_pperp').shape) != 1 or g('_pperp').size != c.shape[2]:
                    raise TransportException("{}: Invalid dimensions of pperp vector. Expected {} elements.".format(coeff, c.shape[2]))
                elif not np.isarray(g('_ppar')) or len(g('_ppar').shape) != 1 or g('_ppar').size != c.shape[3]:
                    raise TransportException("{}: Invalid dimensions of ppar vector. Expected {} elements.".format(coeff, c.shape[3]))
            else:
                raise TransportException("No momentum grid provided for transport coefficient '{}'.".format(coeff))
        else:
            if len(c.shape) != 4:
                raise TransportException("{}: Invalid dimensions of transport coefficient: {}".format(coeff, c.shape))
            elif not np.isarray(g('_t')) or len(g('_t').shape) != 1 or g('_t').size != c.shape[0]:
                raise TransportException("{}: Invalid dimensions of time vector. Expected {} elements.".format(coeff, c.shape[0]))
            elif not np.isarray(g('_r')) or len(g('_r').shape) != 1 or g('_r').size != c.shape[1]:
                raise TransportException("{}: Invalid dimensions of radius vector. Expected {} elements.".format(coeff, c.shape[1]))


class TransportException(DREAMException):
    def __init__(self, msg):
        super().__init__(msg)


