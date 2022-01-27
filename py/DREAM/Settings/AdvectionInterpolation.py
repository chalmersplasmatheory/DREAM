# Routines for setting advection interpolation methods.


AD_INTERP_CENTRED  = 1
AD_INTERP_UPWIND   = 2
AD_INTERP_UPWIND_2ND_ORDER = 3
AD_INTERP_DOWNWIND = 4
AD_INTERP_QUICK    = 5
AD_INTERP_SMART    = 6
AD_INTERP_MUSCL    = 7
AD_INTERP_OSPRE    = 8
AD_INTERP_TCDF     = 9

AD_INTERP_JACOBIAN_LINEAR = 1
AD_INTERP_JACOBIAN_FULL   = 2
AD_INTERP_JACOBIAN_UPWIND = 3


class AdvectionInterpolation:
    

    def __init__(self, kinetic=True,
        ad_int_r=AD_INTERP_CENTRED, ad_int_p1=AD_INTERP_CENTRED,
        ad_int_p2=AD_INTERP_CENTRED, ad_jac_r=AD_INTERP_JACOBIAN_FULL,
        ad_jac_p1=AD_INTERP_JACOBIAN_FULL, ad_jac_p2=AD_INTERP_JACOBIAN_FULL,
        fluxlimiterdamping=1.0):
        """
        Constructor.
        """
        self.kinetic = kinetic

        self.adv_interp_r  = ad_int_r
        self.adv_interp_p1 = ad_int_p1
        self.adv_interp_p2 = ad_int_p2

        self.adv_jac_r  = ad_jac_r
        self.adv_jac_p1 = ad_jac_p1
        self.adv_jac_p2 = ad_jac_p2

        self.fluxlimiterdamping = fluxlimiterdamping


    def setMethod(self, ad_int=None, ad_int_r=AD_INTERP_CENTRED,
        ad_int_p1=AD_INTERP_CENTRED, ad_int_p2=AD_INTERP_CENTRED, ad_jac=None, 
        ad_jac_r=AD_INTERP_JACOBIAN_LINEAR, ad_jac_p1=AD_INTERP_JACOBIAN_LINEAR,
        ad_jac_p2=AD_INTERP_JACOBIAN_LINEAR, fluxlimiterdamping=1.0):
        """
        Sets the interpolation method that is used in the advection terms of
        the equation. To set all three components simultaneously, provide ad_int
        and/or ad_jac. Otherwise the three components can use separate interpolation
        methods.
        
        :param int ad_int:               Interpolation method to use for all coordinates.
        :param int ad_int_r:             Interpolation method to use for the radial coordinate.
        :param int ad_int_p1:            Interpolation method to use for the first momentum coordinate.
        :param int ad_int_p2:            Interpolation method to use for the second momentum coordinate.
        :param int ad_jac:               Jacobian interpolation mode to use for all coordinates.
        :param int ad_jac_r:             Jacobian interpolation mode to use for the radial coordinate.
        :param int ad_jac_p1:            Jacobian interpolation mode to use for the first momentum coordinate.
        :param int ad_jac_p2:            Jacobian interpolation mode to use for the second momentum coordinate.
        :param float fluxlimiterdamping: Damping parameter used to under-relax the interpolation coefficients during non-linear iterations (should be between 0 and 1).
        """
        self.fluxlimiterdamping = fluxlimiterdamping
        if ad_int is not None:
            self.adv_interp_r  = ad_int
            self.adv_interp_p1 = ad_int
            self.adv_interp_p2 = ad_int
        else:
            self.adv_interp_r  = ad_int_r
            self.adv_interp_p1 = ad_int_p1
            self.adv_interp_p2 = ad_int_p2

        if ad_jac is not None:
            self.adv_jac_r  = ad_jac
            self.adv_jac_p1 = ad_jac
            self.adv_jac_p2 = ad_jac
        else:
            self.adv_jac_r  = ad_jac_r
            self.adv_jac_p1 = ad_jac_p1
            self.adv_jac_p2 = ad_jac_p2


    def fromdict(self, data):
        """
        Load settings from the specified dictionary.
        """
        self.adv_interp_r = data['r']
        self.adv_jac_r = data['r_jac']

        if self.kinetic:
            self.adv_interp_p1 = data['p1'][:]
            self.adv_interp_p2 = data['p2'][:]
            self.adv_jac_p1 = data['p1_jac'][:]
            self.adv_jac_p2 = data['p2_jac'][:]

        self.fluxlimiterdamping = data['fluxlimiterdamping'][:]


    def todict(self):
        """
        Convert these settings to a dictionary.
        """
        data = {}

        data['r']  = self.adv_interp_r
        data['r_jac'] = self.adv_jac_r

        if self.kinetic:
            data['p1'] = self.adv_interp_p1
            data['p2'] = self.adv_interp_p2
            data['p1_jac'] = self.adv_jac_p1
            data['p2_jac'] = self.adv_jac_p2

        data['fluxlimiterdamping'] = self.fluxlimiterdamping

        return data


    def verifySettings(self):
        ad_int_opts = [
            AD_INTERP_CENTRED, AD_INTERP_DOWNWIND, AD_INTERP_UPWIND, AD_INTERP_UPWIND_2ND_ORDER, 
            AD_INTERP_QUICK, AD_INTERP_SMART, AD_INTERP_MUSCL, AD_INTERP_OSPRE, AD_INTERP_TCDF
        ]
        if self.adv_interp_r not in ad_int_opts:
            raise EquationException("{}: Invalid radial interpolation coefficient set: {}.".format(self.name, self.adv_interp_r))

        if self.kinetic:
            if self.adv_interp_p1 not in ad_int_opts:
                raise EquationException("{}: Invalid p1 interpolation coefficient set: {}.".format(self.name, self.adv_interp_p1))
            if self.adv_interp_p2 not in ad_int_opts:
                raise EquationException("{}: Invalid p2 interpolation coefficient set: {}.".format(self.name, self.adv_interp_p2))

        if (self.fluxlimiterdamping<0.0) or (self.fluxlimiterdamping>1.0):
            raise EquationException("{}: Invalid flux limiter damping coefficient: {}. Choose between 0 and 1.".format(self.name, self.fluxlimiterdamping))


