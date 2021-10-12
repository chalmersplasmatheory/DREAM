# Settings for the SPI shards (sizes and velocities)

import numpy as np
from scipy.special import kn
from scipy import integrate
from scipy.constants import N_A
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity
import DREAM.Settings.Equations.IonSpecies as Ions



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

ZMolarMassList=[1,1,10]
isotopesMolarMassList=[2,0,0]# 0 means naturally occuring mix
molarMassList=[0.0020141,0.001008,0.020183]# kg/mol

ZSolidDensityList=[1,1,10]
isotopesSolidDensityList=[2,0,0]
solidDensityList=[205.9,86,1444]# kg/m^3

class SPI(UnknownQuantity):
    

    def __init__(self, settings, rp=None, vp=None, xp=None , VpVolNormFactor=1, rclPrescribedConstant=0.01, velocity=VELOCITY_MODE_NONE, ablation=ABLATION_MODE_NEGLECT, deposition=DEPOSITION_MODE_NEGLECT, heatAbsorbtion=HEAT_ABSORBTION_MODE_NEGLECT, cloudRadiusMode=CLOUD_RADIUS_MODE_NEGLECT):
        """
        Constructor.
        
        :param DREAMSettings settings: Parent DREAMSettings object.
        :param numpy.ndarray rp: Initial shard radii.
        :param numpy.ndarray vp: Initial shard velocities (cartesian coordinates)
        :param numpy.ndarray xp: Initial shard positions (cartesian coordinates)
        :param float VpVolNormFactor: Factor used to renormalize the value of VpVol 
                used for calculating the voluma of the flux tubes (eg major radius in 
                case of cylindrical geometry)
        :param float rclPrescribedConstant: Constant, prescribed radius of the neutral 
                cloud surrounding each pellet shard (only applicable if 
                cloudRadiusMode=CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT)
        :param int velocity: Model used for the shard velocities
        :param int ablation: Model used for shard ablation
        :param int deposition: Model used for the deposition of the ablated material
        :param int heatAbsobtion: Model used for absorbtion of heat flowing into the neutral clouds
        :param int cloudRadiusMode: Mode used for calculating the radius of the neutral clouds
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


    def setInitialData(self, rp=None, vp=None, xp=None):

        if rp is not None:
            if np.isscalar(rp):
                self.rp = np.asarray([rp])
            else: self.rp = np.asarray(rp)

        if vp is not None:
            if np.isscalar(vp):
                self.vp = np.asarray([vp])
            else: self.vp = np.asarray(vp)

        if xp is not None:
            if np.isscalar(xp):
                self.xp = np.asarray([xp])
            else: self.xp = np.asarray(xp)
        
    def rpDistrParksStatistical(self,rp,kp):
        """
        Evaluates the shard size distribution function referred to as the 
        'statistical model' in P. Parks 2016 GA report (DOI:10.2172/1344852)
        """
        return kn(0,rp*kp)*kp**2*rp
        
    def sampleRpDistrParksStatistical(self,N,kp):
        """
        Samples N shard radii according to the distribution function 
        given by rpDistrParksStatistical()
        """
        # First we calculate the cdf, and then interpolate the cdf-values 
        # back to the corresponding radii at N randomly chosen points between 0 and 1
        rp_integrate=np.linspace(1e-10/kp,10/kp,5000)
        cdf=integrate.cumtrapz(y=self.rpDistrParksStatistical(rp_integrate,kp),x=rp_integrate)
        return np.interp(np.random.uniform(size=N),np.hstack((0.0,cdf)),rp_integrate)
        
    def setRpParksStatistical(self, nShard, Ninj, Zs, isotopes, molarFractions, ionNames,  opacity_modes = None, add=True, **kwargs):
        """
        sets (or adds) nShard shards with radii distributed accordin to 
        rpDistrParksStatistical(), with the characteristic inverse shard size kp 
        calculated from the given pellet and shattering parameters. Also updates the ion 
        settings with the appropriate molar fractions contributing to each ion species
        
        :param int nShard: Number of shards into which the pellet is shattered
        :param float Ninj: Numbr of particles contained in the pellet
        :param list Zs: List of charge numbers for every ion species the pellet consists of
        :param list isotopes: List of isotopes for every ion species the pellet consists of
        :param numpy.ndarray molarFractions: Molar fraction with which each ion species contribute
        :param list ionNames: List of names for the ion species to be added and connected 
                  to the ablation of this pellet
        :param list opacity_modes: List of opacity modes for every ion species the pellet consists of.
                  If 'None', this argument is omitted when adding the ion species, so that the default 
                  settings (transparent) is used
        :param bool add: If 'True', add the new pellet shards to the existing ones, otherwise 
             existing shards are cleared
             
        :return: the inverse characteristic shard size kp
        """
        
        
        # Calculate solid particle density of the pellet (needed to calculate the 
        # inverse characteristic shard size)
        molarVolume=0
        for iZ in range(len(Zs)):
            for iList in range(len(solidDensityList)):
                if Zs[iZ]==ZSolidDensityList[iList] and isotopes[iZ]==isotopesSolidDensityList[iList]:
                    solidDensityIZ=solidDensityList[iList]
                if Zs[iZ]==ZMolarMassList[iList] and isotopes[iZ]==isotopesMolarMassList[iList]:
                    molarMassIZ=molarMassList[iList]
            
            molarVolume+=molarFractions[iZ]*molarMassIZ/solidDensityIZ
            
        solidParticleDensity=N_A/molarVolume
       
       
        # Calculate inverse characteristic shard size
        kp=(6*np.pi**2*solidParticleDensity*nShard/Ninj)**(1/3)
        
        # Sample the shard sizes and rescale to get exactly the 
        # specified number of particles in the pellet
        rp_init=self.sampleRpDistrParksStatistical(nShard,kp)
        Ninj_obtained=np.sum(4*np.pi*rp_init**(3)/3/molarVolume*N_A)
        rp_init*=(Ninj/Ninj_obtained)**(1/3)       
       
        if add and self.rp is not None:
            self.rp=np.concatenate((self.rp,rp_init))
        else:
            self.rp=rp_init
            
        # Add zeros to the end of SPIMolarFraction for all ion species previously connected to a pellet
        for ion in self.settings.eqsys.n_i.ions:
            SPIMolarFractionPrevious=ion.getSPIMolarFraction()
            if SPIMolarFractionPrevious[0]!=-1:
                ion.setSPIMolarFraction(np.concatenate((SPIMolarFractionPrevious,np.zeros(nShard))))
                
                
        # Add an ion species connected to this pellet to the ion settings
        for iZ in range(len(Zs)):
            
            # SPIMolarFraction must have the smae length as all pellet shard, 
            # not only the pellet which is initiated here, so set the molar fraction 
            # to zero for previously set shards
            SPIMolarFraction=np.zeros(len(self.rp))
            SPIMolarFraction[-nShard:]=molarFractions[iZ]*np.ones(nShard)
            if opacity_modes is not None:
                self.settings.eqsys.n_i.addIon(name=ionNames[iZ], n=1e0, Z=Zs[iZ], isotope=isotopes[iZ], opacity_mode=opacity_modes[iZ], iontype=Ions.IONS_DYNAMIC_NEUTRAL, SPIMolarFraction=SPIMolarFraction,**kwargs)
            
            else:
                self.settings.eqsys.n_i.addIon(name=ionNames[iZ], n=1e0, Z=Zs[iZ], isotope=isotopes[iZ], iontype=Ions.IONS_DYNAMIC_NEUTRAL, SPIMolarFraction=SPIMolarFraction,**kwargs)
            
        return kp
        
    def setShardPositionSinglePoint(self, nShard,shatterPoint,add=True):
        """
        Sets self.xp to a vector of the (x,y,z)-coordinates of nShard initial
        pellet shard positions starting from the single point shatterPoint
        
        :param int nShard: Number of shards 
        :param numpy.ndarray shatterPoint: (x,y,z)-coordinates for the starting point of the shards to be set
        :param bool add: If 'True', add the new pellet shard positions to the existing ones, otherwise 
             existing shards are cleared
        """
        if add and self.xp is not None:
            self.xp=np.concatenate((self.xp,np.tile(shatterPoint,nShard)))
        else:
            self.xp=np.tile(shatterPoint,nShard)
            
    def setShardVelocitiesUniform(self, nShard,abs_vp_mean,abs_vp_diff,alpha_max,nDim=2,add=True, shards=None):
        """
        Sets self.vp to a vector storing the (x,y,z)-components of nShard shard velosities,
        assuming a uniform velocity distribution over a nDim-dimensional cone whose axis
        is anti-parallell to the x-axis. TODO: implement support for an arbitrary axis?
        
        :param int nShard: Number of shards
        :param float abs_vp_mean: Mean of the magnitude of the shard velocities
        :param float abs_vp_diff: width of the uniform distribution of the magnitude of the shard velocities
        :param float alpha_max: Span of divergence angle (ie twice the opening angle of the cone)
        :param int nDim: number of dimensions into which the shards should be spread
        :param bool add: If 'True', add the new pellet shard velocities to the existing ones, otherwise 
             existing shards are cleared
        :param slice shards: indices of existing shards whose velocities should be updated. If not 'None', 
                add is set to 'False' and nShard is set to the number of indices to be updated
        """
        
        if shards is not None:
        	nShard=len(self.vp[shards])
        	add=False
        
        # Sample magnitude of velocities
        abs_vp_init=(abs_vp_mean+abs_vp_diff*(-1+2*np.random.uniform(size=nShard)))
        
        # Sample directions uniformly over a nDim-dimensional cone and set the velocity vectors
        vp_init=np.zeros(3*nShard)
        if nDim==1:
            # in 1D, the "cone" simply becomes a straight line
            vp_init[0::3]=-abs_vp_init
            
        elif nDim==2:
            # in 2D, the cone becomes a circle sector
            alpha=alpha_max*(-1+2*np.random.uniform(size=nShard))
            vp_init[0::3]=-abs_vp_init*np.cos(alpha)
            vp_init[1::3]=abs_vp_init*np.sin(alpha)
            
        elif nDim==3:
            # The solid angle covered by the part of the cone between alpa and d(alpha) 
            # is proportional to sin(alpha), and the normalised probability distribution 
            # becomes f(alpha)=sin(alpha)/(1-cos(alpha_max/2)). We sample from this
            # distribution by applying the inverse cdf to uniformly drawn numbers
            # between 0 and 1
            alpha=np.arccos(1-np.random.uniform(size=nShard)*(1-np.cos(alpha_max/2)))
            
            # The angle in the yz-plane is simply drawn randomly
            phi=2*np.pi*np.random.uniform(size=nShard)
            
            # Finally calculate the velocity vectors
            vp_init[0::3]=-abs_vp_init*np.cos(alpha)
            vp_init[1::3]=abs_vp_init*np.sin(alpha)*np.cos(phi)
            vp_init[2::3]=abs_vp_init*np.sin(alpha)*np.sin(phi)
            
        else:
            raise EquationException("spi: Invalid number of dimensions into which the pellet shards are spread")
            
        if add and self.vp is not None:
            self.vp=np.concatenate((self.vp,vp_init))
        elif shards is not None:
        	# Pick out the components of the stored shard velocities...
        	vpx=self.vp[0::3]
        	vpy=self.vp[1::3]
        	vpz=self.vp[2::3]
        	
        	# ... Change the velocities of the shards specified in the input...
        	vpx[shards]=vp_init[0::3]
        	vpy[shards]=vp_init[1::3]
        	vpz[shards]=vp_init[2::3]
        	
        	# ...and finallyset the stored velocities to the updated ones
        	self.vp[0::3]=vpx
        	self.vp[1::3]=vpy
        	self.vp[2::3]=vpz
        else:
            self.vp=vp_init
            
    def setParamsVallhagenMSc(self, nShard, Ninj, Zs, isotopes, molarFractions, ionNames, shatterPoint, abs_vp_mean,abs_vp_diff,alpha_max,nDim=2, add=True, opacity_modes = None, **kwargs):
        """
        Wrapper for setRpParksStatistical(), setShardPositionSinglePoint() and setShardVelocitiesUniform(),
        which combined are used to set up an SPI-scenario similar to those in Oskar Vallhagens MSc thesis
        (available at https://hdl.handle.net/20.500.12380/302296)
        """
        
        kp=self.setRpParksStatistical(nShard, Ninj, Zs, isotopes, molarFractions, ionNames, opacity_modes, add, **kwargs)
        self.setShardPositionSinglePoint(nShard,shatterPoint,add)
        self.setShardVelocitiesUniform(nShard,abs_vp_mean,abs_vp_diff,alpha_max,nDim,add)
        return kp
        
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
        
        # If no SPI settings have been given, set everything to zero (to avoid a DREAMIOException)
        # Before this stage it is usefull to use None to indicate if any SPI settings have been made yet,
        # to know if there are any previous shards to add the new ones to, so therefore
        # we don't set this default setting until this stage
        if self.rp is None:
            self.rp=np.array([0])
        if self.vp is None:
            self.vp=np.array([0,0,0])
        if self.xp is None:
            self.xp=np.array([0,0,0])
            
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
            raise EquationException("spi: Invalid value assigned to 'velocity'. Expected integer.")
        if type(self.ablation) != int:
            raise EquationException("spi: Invalid value assigned to 'ablation'. Expected integer.")
        if type(self.deposition) != int:
            raise EquationException("spi: Invalid value assigned to 'deposition'. Expected integer.")
        if type(self.heatAbsorbtion) != int:
            raise EquationException("spi: Invalid value assigned to 'heatAbsorbtion'. Expected integer.")



    def verifySettingsPrescribedInitialData(self):
        if vp.size!=3*rp.size:
            raise EquationException("Missmatch in size of initial data arrays for rp and vp. Expected vp to have a size 3 times the size of rp")
        if xp.size!=3*rp.size:
            raise EquationException("Missmatch in size of initial data arrays for rp and xp. Expected xp to have a size 3 times the size of rp")
        
