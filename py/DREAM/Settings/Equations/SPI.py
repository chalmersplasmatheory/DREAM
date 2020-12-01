# Settings for the runaway electron density

import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity



VELOCITY_MODE_NONE=1
VELOCITY_MODE_PRESCRIBED=2

ABLATION_MODE_NEGLECT=1
ABLATION_MODE_FLUID_NGS=2
ABLATION_MODE_KINETIC_NGS=3

DEPOSITION_MODE_NEGLECT=1
DEPOSITION_MODE_LOCAL=2
DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE=3
DEPOSITION_MODE_LOCAL_GAUSSIAN=4

HEAT_ABSORBTION_MODE_NEGLECT=1
HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS=2
HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN=3

CLOUD_RADIUS_MODE_NEGLECT=1
CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT=2
CLOUD_RADIUS_MODE_SELFCONSISTENT=3

class SPI(UnknownQuantity):
    

    def __init__(self, settings, rp=0, vp=0, xp=0 , VpVolNormFactor=1, rclPrescribedConstant=0.01, velocity=VELOCITY_MODE_NONE, ablation=ABLATION_MODE_NEGLECT, deposition=DEPOSITION_MODE_NEGLECT, heatAbsorbtion=HEAT_ABSORBTION_MODE_NEGLECT, cloudRadiusMode=CLOUD_RADIUS_MODE_NEGLECT):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.velocity           = velocity
        self.ablation           = ablation
        self.deposition         = deposition
        self.heatAbsorbtion     = heatAbsorbtion
        self.cloudRadiusMode    = cloudRadiusMode
        self.VpVolNormFactor    = VpVolNormFactor
        self.rclPrescribedConstant = rclPrescribedConstant

        self.rp       = None
        self.vp       = None
        self.xp       = None
        self.setInitialData(rp=rp, vp=vp, xp=xp)


    def setInitialData(self, rp, vp, xp):

        if np.isscalar(rp):
            self.rp = np.asarray([rp])
        else: self.rp = np.asarray(rp)

        if np.isscalar(vp):
            self.vp = np.asarray([vp])
        else: self.vp = np.asarray(vp)

        if np.isscalar(xp):
            self.xp = np.asarray([xp])
        else: self.xp = np.asarray(xp)

    def setVpVolNormFactor(self,VpVolNormFactor):
        self.VpVolNormFactor=VpVolNormFactor

    def setRclPrescribedConstant(self,rclPrescribedConstant):
        self.rclPrescribedConstant=rclPrescribedConstant


    def setVelocity(self, velocity):
        """
        Specifies mode to calculate shard velocities.
        """
        self.velocity = int(velocity)


    def setAblation(self, ablation):
        """
        Specifies which model to use for calculating the
        ablation rate.
        """
        self.ablation = int(ablation)

    def setDeposition(self, deposition):
        """
        Specifies which model to use for calculating the
        deposition of ablated material.
        """
        self.deposition = int(deposition)

    def setHeatAbsorbtion(self, heatAbsorbtion):
        """
        Specifies which model to use for calculating the
        heat absorbtion in the neutral pellet cloud
        """
        self.heatAbsorbtion = int(heatAbsorbtion)

    def setCloudRadiusMode(self, cloudRadiusMode):
        """
        Specifies which model to use for calculating the
        radius of the the neutral pellet cloud
        """
        self.cloudRadiusMode = int(cloudRadiusMode)


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        self.velocity       = data['velocity']
        self.ablation       = data['ablation']
        self.deposition     = data['deposition']
        self.heatAbsorbtion = data['heatAbsorbtion']
        self.cloudRadiusMode = data['cloudRadiusMode']

        self.VpVolNormFactor = data['VpVolNormFactor']
        self.rclPrescribedConstant = data['rclPrescribedConstant']
        self.rp              = data['init']['rp']
        self.vp              = data['init']['vp']
        self.xp              = data['init']['xp']


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this SPI object.
        """
        data = {
            'velocity': self.velocity,
            'ablation': self.ablation,
            'deposition': self.deposition,
            'heatAbsorbtion': self.heatAbsorbtion,
            'cloudRadiusMode': self.cloudRadiusMode,
            'VpVolNormFactor': self.VpVolNormFactor,
            'rclPrescribedConstant': self.rclPrescribedConstant
        }
        data['init'] = {
                'rp': self.rp,
                'vp': self.vp,
                'xp': self.xp
        }

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if type(self.velocity) != int:
            raise EquationException("n_re: Invalid value assigned to 'velocity'. Expected integer.")
        if type(self.ablation) != int:
            raise EquationException("n_re: Invalid value assigned to 'ablation'. Expected integer.")
        if type(self.deposition) != int:
            raise EquationException("n_re: Invalid value assigned to 'deposition'. Expected integer.")
        if type(self.heatAbsorbtion) != int:
            raise EquationException("n_re: Invalid value assigned to 'heatAbsorbtion'. Expected integer.")



    def verifySettingsPrescribedInitialData(self):
        if vp.size!=3*rp.size:
            raise EquationException("Missmatch in size of initial data arrays for rp and vp. Expected vp to have a size 3 times the size of rp")
        if xp.size!=3*rp.size:
            raise EquationException("Missmatch in size of initial data arrays for rp and xp. Expected xp to have a size 3 times the size of rp")
        
