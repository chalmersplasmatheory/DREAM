# General transport settings which can be added to any UnknownQuantity


import numpy as np
from .. DREAMException import DREAMException


TRANSPORT_NONE = 1
TRANSPORT_PRESCRIBED = 2
TRANSPORT_RECHESTER_ROSENBLUTH = 3
TRANSPORT_SVENSSON = 4
TRANSPORT_FROZEN_CURRENT = 5

INTERP3D_NEAREST = 0
INTERP3D_LINEAR  = 1

INTERP1D_NEAREST = 0
INTERP1D_LINEAR  = 1

SVENSSON_INTERP1D_PARAM_TIME = 1
SVENSSON_INTERP1D_PARAM_IP   = 2

FROZEN_CURRENT_MODE_DISABLED = 1
FROZEN_CURRENT_MODE_CONSTANT = 2
FROZEN_CURRENT_MODE_BETAPAR = 3

BC_CONSERVATIVE = 1     # Assume no flux through r=rmax
BC_F_0 = 2              # Assume f=0 outside the plasma
BC_DF_CONST = 3         # Assume that df/dr is constant on the plasma boundary


class TransportSettings:
    

    def __init__(self, kinetic=False):
        """
        Constructor.

        :param bool kinetic: If ``True``, the coefficient will be assumed kinetic (4D). Otherwise fluid (2D).
        """
        self.kinetic = kinetic
        self.type    = TRANSPORT_NONE


        # Prescribed advection
        self.ar           = None
        self.ar_t         = None
        self.ar_r         = None
        self.ar_p         = None
        self.ar_xi        = None
        self.ar_ppar      = None
        self.ar_pperp     = None
        self.ar_interp3d  = None
        
        # Prescribed diffusion
        self.drr          = None
        self.drr_t        = None
        self.drr_r        = None
        self.drr_p        = None
        self.drr_xi       = None
        self.drr_ppar     = None
        self.drr_pperp    = None
        self.drr_interp3d = None

        # Svensson pstar
        self.pstar          = None
        self.interp1d_param = SVENSSON_INTERP1D_PARAM_TIME 
        
        # Svensson advection
        self.s_ar           = None
        self.s_ar_r         = None
        self.s_ar_t         = None
        self.s_ar_p         = None
        self.s_ar_xi        = None
        self.s_ar_ppar      = None
        self.s_ar_pperp     = None
        self.s_ar_interp3d  = None
        self.s_ar_interp1d  = None
        
        # Svensson diffusion
        self.s_drr          = None
        self.s_drr_r        = None
        self.s_drr_t        = None
        self.s_drr_p        = None
        self.s_drr_xi       = None
        self.s_drr_ppar     = None
        self.s_drr_pperp    = None
        self.s_drr_interp3d = None
        self.s_drr_interp1d = None

        # Rechester-Rosenbluth (diffusive) transport
        self.dBB   = None
        self.dBB_t = None
        self.dBB_r = None

        # Frozen current mode transport
        self.frozen_current_mode = FROZEN_CURRENT_MODE_DISABLED
        self.frozen_current_Ip_presc = None
        self.frozen_current_Ip_presc_t = None
        self.frozen_current_D_I_min = 0
        self.frozen_current_D_I_max = 1000

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
        self._prescribeCoefficient('drr', coeff=drr, t=t, r=r, p=p, xi=xi, ppar=ppar, pperp=pperp)


    def setSvenssonPstar(self,pstar):
        """
        Set the lower momentum bound for the runaway, radial transport, region.
        """
        self.pstar=float(pstar)
    

    def setSvenssonInterp1dParam(self, interp1d_param=SVENSSON_INTERP1D_PARAM_TIME):
        """
        Set the lower momentum bound for the runaway, radial transport, region.
        """
        self.interp1d_param = int(interp1d_param)


    def setBoundaryCondition(self, bc=None):
        """
        Set the type of boundary condition. (Default is BC_CONSERVATIVE)
        """
        self.boundarycondition = bc
    

    def setSvenssonAdvection(self, ar, t=None, Ip=None, r=None, p=None, xi=None, ppar=None, pperp=None, interp3d=INTERP3D_LINEAR, interp1d=INTERP1D_LINEAR):
        r"""
        Set the Svensson advection coefficient to use.

        :param ar:       Advection coefficient, :math:`A_r(t,r,\xi_0,p)` or :math:`A_r(I_p,r,\xi_0,p)`.
        :param t:        Time vector for which ``ar`` is defined (if ``Ip`` is not provided).
        :param Ip:       Plasma current vector for which ``ar`` is defined (if ``t`` is not provided).
        :param r:        Radial grid vector for which ``ar`` is defined.
        :param p:        Momentum grid vector for which ``ar`` is defined.
        :param xi:       Pitch grid vector for which ``ar`` is defined.
        :param interp3d: Interpolation method to use when interpolating in (r,xi,p) part of coefficient.
        :param interp1d: Interpolation method to use when interpolating in time/Ip variable.
        """
        if self.interp1d_param == SVENSSON_INTERP1D_PARAM_TIME:
            if t is not None:
                self._prescribeCoefficient('s_ar', coeff=ar, t=t, r=r, p=p, xi=xi, ppar=ppar, pperp=pperp,interp3d=interp3d,override_kinetic=True)
            else: 
                raise TransportException('interp1d_param has been set to "time", but no time variable was given.')
        elif self.interp1d_param == SVENSSON_INTERP1D_PARAM_IP:
            if Ip is not None:
                self._prescribeCoefficient('s_ar', coeff=ar, t=Ip, r=r, p=p, xi=xi, ppar=ppar, pperp=pperp,interp3d=interp3d,override_kinetic=True)
            else:
                raise TransportException('interp1d_param has been set to "Ip", but no plasma-current variable was given.')
        else:
            raise TransportException('interp1d_param has not been set or is invalid. It must be set before setting the Svensson transport coefficients.')
        self.type = TRANSPORT_SVENSSON
        self.s_ar_interp1d = interp1d
    

    def setSvenssonDiffusion(self, drr, t=None, Ip=None, r=None, p=None, xi=None, ppar=None, pperp=None,interp3d=INTERP3D_LINEAR, interp1d=INTERP1D_LINEAR):
        r"""
        Set the Svensson diffusion coefficient to use.

        :param drr:      Diffusion coefficient, :math:`D_{rr}(t,r,\xi_0,p)` or :math:`D_{rr}(I_p,r,\xi_0,p)`.
        :param t:        Time vector for which ``drr`` is defined (if ``Ip`` is not provided).
        :param Ip:       Plasma current vector for which ``drr`` is defined (if ``t`` is not provided).
        :param r:        Radial grid vector for which ``drr`` is defined.
        :param p:        Momentum grid vector for which ``drr`` is defined.
        :param xi:       Pitch grid vector for which ``drr`` is defined.
        :param interp3d: Interpolation method to use when interpolating in (r,xi,p) part of coefficient.
        :param interp1d: Interpolation method to use when interpolating in time/Ip variable.
        """
        if self.interp1d_param == SVENSSON_INTERP1D_PARAM_TIME:
            if t is not None:
                self._prescribeCoefficient('s_drr', coeff=drr, t=t, r=r, p=p, xi=xi, ppar=ppar, pperp=pperp,interp3d=interp3d,override_kinetic=True)
            else: 
                raise TransportException('interp1d_param has been set to "time", but no time variable was given.')
        elif self.interp1d_param == SVENSSON_INTERP1D_PARAM_IP:
            if Ip is not None:
                self._prescribeCoefficient('s_drr', coeff=drr, t=Ip, r=r, p=p, xi=xi, ppar=ppar, pperp=pperp,interp3d=interp3d,override_kinetic=True)
            else:
                raise TransportException('interp1d_param has been set to "Ip", but no plasma-current variable was given.')
        else:
            raise TransportException('interp1d_param has not been set or is invalid. It must be set before setting the Svensson transport coefficients.')
        self.type = TRANSPORT_SVENSSON
        self.s_drr_interp1d = interp1d


    def _prescribeCoefficient(self, name, coeff, t=None, r=None, p=None, xi=None, ppar=None, pperp=None,interp3d=INTERP3D_LINEAR, override_kinetic=False):
        """
        General method for prescribing an advection or diffusion coefficient.
        """
        self.type = TRANSPORT_PRESCRIBED

        setattr(self, name+'_interp3d', interp3d)

        if np.isscalar(coeff):
            r = np.array([0])
            t = np.array([0])
            p = np.array([0])
            xi = np.array([0])

            if self.kinetic or override_kinetic:
                coeff = coeff * np.ones((1,)*4)
            else:
                coeff = coeff * np.ones((1,)*2)

        r = np.asarray(r)
        t = np.asarray(t)
        
        if r.ndim != 1: r = np.reshape(r, (r.size,))
        if t.ndim != 1: t = np.reshape(t, (t.size,))

        if (self.kinetic == False and not override_kinetic) and len(coeff.shape) == 2:
            setattr(self, name, coeff)
            setattr(self, name+'_r', r)
            setattr(self, name+'_t', t)
        elif (self.kinetic == True or override_kinetic) and len(coeff.shape) == 4:
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
                setattr(self, name+'_p', p)
                setattr(self, name+'_xi', xi)
            else:
                setattr(self, name+'_ppar', ppar)
                setattr(self, name+'_pperp', pperp)
        else:
            raise TransportException("Invalid dimensions of prescribed coefficient: {}. Expected {} dimensions.".format(coeff.shape, 4 if (self.kinetic or override_kinetic) else 2))

            
    def setMagneticPerturbation(self, dBB, t=None, r=None):
        """
        Prescribes the evolution of the magnetic perturbation level (dB/B).

        :param dBB: Magnetic perturbation level.
        :param t:   Time grid on which the perturbation is defined.
        :param r:   Radial grid on which the perturbation is defined.
        """
        self.type = TRANSPORT_RECHESTER_ROSENBLUTH

        if np.isscalar(dBB):
            dBB = dBB * np.ones((1,1))
            r = np.array([0])
            t = np.array([0])

        r = np.asarray(r)
        t = np.asarray(t)

        if r.ndim != 1: r = np.reshape(r, (r.size,))
        if t.ndim != 1: t = np.reshape(t, (t.size,))

        self.dBB_r = r
        self.dBB_t = t
        self.dBB   = dBB


    def setFrozenCurrentMode(self, mode, Ip_presc, Ip_presc_t=0, D_I_min=0, D_I_max=1000):
        """
        Enable the frozen current mode and specify the target plasma current.
        """
        self.type = TRANSPORT_FROZEN_CURRENT
        self.frozen_current_mode = mode

        if np.isscalar(Ip_presc) or Ip_presc.ndim == 0:
            Ip_presc = np.array([Ip_presc])
            Ip_presc_t = np.array([0])
        if not isinstance(Ip_presc_t, np.ndarray):
            Ip_presc_t = np.array([Ip_presc_t])

        if Ip_presc_t.ndim != 1:
            Ip_presc_t = np.reshape(Ip_presc_t, (Ip_presc_t.size,))
            
        self.frozen_current_Ip_presc = Ip_presc
        self.frozen_current_Ip_presc_t = Ip_presc_t
        self.frozen_current_D_I_min = D_I_min
        self.frozen_current_D_I_max = D_I_max


    def setBoundaryCondition(self, bc):
        """
        Set the boundary condition to use for the transport.
        """
        self.boundarycondition = bc


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        # Prescribed advection
        self.ar = None
        self.ar_r = None
        self.ar_t = None
        self.ar_p = None
        self.ar_xi = None
        self.ar_ppar = None
        self.ar_pperp = None
        self.ar_interp3d =None
        
        # Prescribed diffusion
        self.drr = None
        self.drr_r = None
        self.drr_t = None
        self.drr_p = None
        self.drr_xi = None
        self.drr_ppar = None
        self.drr_pperp = None
        self.drr_interp3d =None
        
        # Svensson pstar
        self.pstar          = None
        self.interp1d_param = None
        
        # Svensson advection
        self.s_ar           = None
        self.s_ar_r         = None
        self.s_ar_t         = None
        self.s_ar_p         = None
        self.s_ar_xi        = None
        self.s_ar_ppar      = None
        self.s_ar_pperp     = None
        self.s_ar_interp3d  = None
        self.s_ar_interp1d  = None
        

        # Svensson diffusion
        self.s_drr          = None
        self.s_drr_r        = None
        self.s_drr_t        = None
        self.s_drr_p        = None
        self.s_drr_xi       = None
        self.s_drr_ppar     = None
        self.s_drr_pperp    = None
        self.s_drr_interp3d = None
        self.s_drr_interp1d = None


        # Rechester--Rosenbluth
        self.dBB = None
        self.dBB_r = None
        self.dBB_t = None

        if 'type' in data:
            self.type = data['type']

        if 'boundarycondition' in data:
            self.boundarycondition = data['boundarycondition']

        if 'ar' in data:
            self.ar = data['ar']['x']
            self.ar_r  = data['ar']['r']
            self.ar_t  = data['ar']['t']

            if 'interp3d' in data['ar']:
                self.ar_interp3d = data['ar']['interp3d']

            if self.kinetic:
                if 'p' in data['ar']: self.ar_p = data['ar']['p']
                if 'xi' in data['ar']: self.ar_xi = data['ar']['xi']
                if 'ppar' in data['ar']: self.ar_ppar = data['ar']['ppar']
                if 'pperp' in data['ar']: self.ar_pperp = data['ar']['pperp']

        if 'drr' in data:
            self.drr = data['drr']['x']
            self.drr_r  = data['drr']['r']
            self.drr_t  = data['drr']['t']

            if 'interp3d' in data['drr']:
                self.drr_interp3d = data['drr']['interp3d']

            if self.kinetic:
                if 'p' in data['drr']: self.drr_p = data['drr']['p']
                if 'xi' in data['drr']: self.drr_xi = data['drr']['xi']
                if 'ppar' in data['drr']: self.drr_ppar = data['drr']['ppar']
                if 'pperp' in data['drr']: self.drr_pperp = data['drr']['pperp']

        if 'pstar' in data:
            self.pstar = float(data['pstar'])
            
        if 'interp1d_param' in data:
            self.interp1d_param = data['interp1d_param']
            
        if 's_ar' in data:
            self.s_ar                 = data['s_ar']['x']
            self.s_ar_r               = data['s_ar']['r']
            self.s_ar_t               = data['s_ar']['t']
            self.s_ar_interp3d        = data['s_ar']['interp3d']
            self.s_ar_interp1d        = data['s_ar']['interp1d']
            
            if 'p' in data['s_ar']:     self.s_ar_p     = data['s_ar']['p']
            if 'xi' in data['s_ar']:    self.s_ar_xi    = data['s_ar']['xi']
            if 'ppar' in data['s_ar']:  self.s_ar_ppar  = data['s_ar']['ppar']
            if 'pperp' in data['s_ar']: self.s_ar_pperp = data['s_ar']['pperp']

        if 's_drr' in data:
            self.s_drr                 = data['s_drr']['x']
            self.s_drr_r               = data['s_drr']['r']
            self.s_drr_t               = data['s_drr']['t']
            self.s_drr_interp3d        = data['s_drr']['interp3d']
            self.s_drr_interp1d        = data['s_drr']['interp1d']

            if 'p' in data['s_drr']:     self.s_drr_p     = data['s_drr']['p']
            if 'xi' in data['s_drr']:    self.s_drr_xi    = data['s_drr']['xi']
            if 'ppar' in data['s_drr']:  self.s_drr_ppar  = data['s_drr']['ppar']
            if 'pperp' in data['s_drr']: self.s_drr_pperp = data['s_drr']['pperp']

        if 'dBB' in data:
            self.dBB   = data['dBB']['x']
            self.dBB_r = data['dBB']['r']
            self.dBB_t = data['dBB']['t']

        if 'frozen_current_mode' in data:
            self.frozen_current_mode = int(data['frozen_current_mode'])
        if 'D_I_min' in data:
            self.frozen_current_D_I_min = float(data['D_I_min'])
        if 'D_I_max' in data:
            self.frozen_current_D_I_max = float(data['D_I_max'])
        if 'I_p_presc' in data:
            self.frozen_current_Ip_presc = data['I_p_presc']['x']
            self.frozen_current_Ip_presc_t = data['I_p_presc']['t']


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
                't': self.ar_t,
                'interp3d': self.ar_interp3d
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
                't': self.drr_t,
                'interp3d': self.drr_interp3d
            }

            if self.kinetic:
                if self.drr_p is not None:
                    data['drr']['p'] = self.drr_p
                    data['drr']['xi'] = self.drr_xi
                else:
                    data['drr']['ppar'] = self.drr_ppar
                    data['drr']['pperp'] = self.drr_pperp

        
        # Svensson pstar
        if self.type == TRANSPORT_SVENSSON and self.pstar is not None:
            data['pstar'] = self.pstar

        # Svensson 1d interpolatino method
        if self.type == TRANSPORT_SVENSSON and self.interp1d_param is not None:
            data['interp1d_param'] =  self.interp1d_param
        
        # Svensson Advection?
        if self.type == TRANSPORT_SVENSSON and self.s_ar is not None:
            data['s_ar'] = {
                'x': self.s_ar,
                'r': self.s_ar_r,
                't': self.s_ar_t,
                'interp3d': self.s_ar_interp3d,
                'interp1d': self.s_ar_interp1d,
            }

            if self.s_ar_p is not None:
                data['s_ar']['p'] = self.s_ar_p
                data['s_ar']['xi'] = self.s_ar_xi
            else:
                data['s_ar']['ppar'] = self.s_ar_ppar
                data['s_ar']['pperp'] = self.s_ar_pperp

        # Svensson Diffusion?
        if self.type == TRANSPORT_SVENSSON and self.s_drr is not None:
            data['s_drr'] = {
                'x': self.s_drr,
                'r': self.s_drr_r,
                't': self.s_drr_t,
                'interp3d': self.s_drr_interp3d,
                'interp1d': self.s_drr_interp1d,
            }

            if self.s_drr_p is not None:
                data['s_drr']['p'] = self.s_drr_p
                data['s_drr']['xi'] = self.s_drr_xi
            else:
                data['s_drr']['ppar'] = self.s_drr_ppar
                data['s_drr']['pperp'] = self.s_drr_pperp

        
        if self.type == TRANSPORT_RECHESTER_ROSENBLUTH and self.dBB is not None:
            data['dBB'] = {
                'x': self.dBB,
                'r': self.dBB_r,
                't': self.dBB_t
            }

        data['frozen_current_mode'] = self.frozen_current_mode
        data['D_I_min'] = self.frozen_current_D_I_min
        data['D_I_max'] = self.frozen_current_D_I_max
        if self.frozen_current_Ip_presc is not None:
            data['I_p_presc'] = {
                'x': self.frozen_current_Ip_presc,
                't': self.frozen_current_Ip_presc_t
            }

        return data


    def verifySettings(self):
        """
        Verify that the settings are consistent.
        """
        if self.type == TRANSPORT_NONE:
            pass
        elif self.type == TRANSPORT_PRESCRIBED:
            self.verifySettingsCoefficient('ar')
            self.verifySettingsCoefficient('drr')
            self.verifyBoundaryCondition()
        elif self.type == TRANSPORT_SVENSSON:
            self.verifySettingsCoefficient('s_ar',override_kinetic=True)
            self.verifySettingsCoefficient('s_drr',override_kinetic=True)
            if self.pstar is None or type(self.pstar) != float:
                raise TransportException("pstar not defined or wrong type.")
            elif self.pstar<=0:
                raise TransportException("pstar = %0.3f <= 0 not allowed." % self.pstar)
            
            self.verifyBoundaryCondition() 
        elif self.type == TRANSPORT_RECHESTER_ROSENBLUTH:
            self.verifySettingsRechesterRosenbluth()
            self.verifyBoundaryCondition()
        elif self.type == TRANSPORT_FROZEN_CURRENT:
            self.verifyFrozenCurrent()
            self.verifyBoundaryCondition()
        else:
            raise TransportException("Unrecognized transport type: {}".format(self.type))


    def verifyBoundaryCondition(self):
        """
        Verify that the boundary condition has been correctly configured.
        """
        bcs = [BC_CONSERVATIVE, BC_F_0, BC_DF_CONST]
        if self.boundarycondition not in bcs:
            raise TransportException("Invalid boundary condition specified for transport: {}".format(self.boundarycondition))


    def verifySettingsCoefficient(self, coeff, override_kinetic=False):
        """
        Verify consistency of the named prescribed transport coefficient.
        """
        g = lambda v : self.__dict__[coeff+v]
        c = g('')

        if c is None: return

        if self.kinetic or override_kinetic:
            if c.ndim != 4:
                raise TransportException("{}: Invalid dimensions of transport coefficient: {}".format(coeff, c.shape))
            elif g('_t').ndim != 1 or g('_t').size != c.shape[0]:
                raise TransportException("{}: Invalid dimensions of time vector. Expected {} elements.".format(coeff, c.shape[0]))
            elif g('_r').ndim != 1 or g('_r').size != c.shape[1]:
                raise TransportException("{}: Invalid dimensions of radius vector. Expected {} elements.".format(coeff, c.shape[1]))

            if g('_interp3d') not in [INTERP3D_LINEAR, INTERP3D_NEAREST]:
                raise TransportException("{}: Invalid value assigned to interp3d.".format(coeff))

            if coeff+'v' in self.__dict__:
                if g('_interp1d') not in [INTERP1D_LINEAR, INTERP1D_NEAREST]:
                    raise TransportException("{}: Invalid value assigned to interp1d.".format(coeff))

            if g('_p') is not None or g('_xi') is not None:
                if g('_xi').ndim != 1 or g('_xi').size != c.shape[2]:
                    raise TransportException("{}: Invalid dimensions of xi vector. Expected {} elements.".format(coeff, c.shape[2]))
                elif g('_p').ndim != 1 or g('_p').size != c.shape[3]:
                    raise TransportException("{}: Invalid dimensions of p vector. Expected {} elements.".format(coeff, c.shape[3]))
            elif g('_ppar') is not None or g('_pperp') is not None:
                if g('_pperp').ndim != 1 or g('_pperp').size != c.shape[2]:
                    raise TransportException("{}: Invalid dimensions of pperp vector. Expected {} elements.".format(coeff, c.shape[2]))
                elif g('_ppar').ndim != 1 or g('_ppar').size != c.shape[3]:
                    raise TransportException("{}: Invalid dimensions of ppar vector. Expected {} elements.".format(coeff, c.shape[3]))
            else:
                raise TransportException("No momentum grid provided for transport coefficient '{}'.".format(coeff))
        else:
            if c.ndim != 2:
                raise TransportException("{}: Invalid dimensions of transport coefficient: {}".format(coeff, c.shape))
            elif g('_t').ndim != 1 or g('_t').size != c.shape[0]:
                raise TransportException("{}: Invalid dimensions of time vector. Expected {} elements.".format(coeff, c.shape[0]))
            elif g('_r').ndim != 1 or g('_r').size != c.shape[1]:
                raise TransportException("{}: Invalid dimensions of radius vector. Expected {} elements.".format(coeff, c.shape[1]))

    def verifySettingsRechesterRosenbluth(self):
        """
        Verify consistency of the Rechester-Rosenbluth transport settings.
        """
        if self.dBB.ndim != 2:
            raise TransportException("Rechester-Rosenbluth: Invalid dimensions of transport coefficient: {}".format(self.dBB.shape))
        elif self.dBB_t.ndim != 1 or self.dBB_t.size != self.dBB.shape[0]:
            raise TransportException("Rechester-Rosenbluth: Invalid dimensions of time vector. Expected {} elements.".format(self.dBB.shape[0]))
        elif self.dBB_r.ndim != 1 or self.dBB_r.size != self.dBB.shape[1]:
            raise TransportException("Rechester-Rosenbluth: Invalid dimensions of radius vector. Expected {} elements.".format(self.dBB.shape[1]))


    def verifyFrozenCurrent(self):
        """
        Verify consistency of the frozen current mode settings.
        """
        if self.frozen_current_mode not in [FROZEN_CURRENT_MODE_DISABLED,FROZEN_CURRENT_MODE_CONSTANT,FROZEN_CURRENT_MODE_BETAPAR]:
            raise TransportException(f"Frozen current mode: Invalid type of transport operator: {self.frozen_current_mode}.")
        elif self.frozen_current_Ip_presc.ndim != 1:
            raise TransportException(f"Frozen current mode: Invalid dimensions of prescribed plasma current: {self.frozen_current_Ip_presc.shape}.")
        elif self.frozen_current_Ip_presc_t.ndim != 1:
            raise TransportException(f"Frozen current mode: Invalid dimensions of prescribed plasma current timebase: {self.frozen_current_Ip_presc_t.shape}.")
        elif self.frozen_current_Ip_presc.shape[0] != self.frozen_current_Ip_presc_t.shape[0]:
            raise TransportException(f"Frozen current mode: Dimension mismatch in prescribed plasma current and its time vector: {self.frozen_current_Ip_presc.shape[0]} =/= {self.frozen_current_Ip_presc_t.shape[0]}.")


class TransportException(DREAMException):
    def __init__(self, msg):
        super().__init__(msg)


