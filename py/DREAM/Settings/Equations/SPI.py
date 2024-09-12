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

ABL_IONIZ_MODE_NEUTRAL = 1
ABL_IONIZ_MODE_SINGLY_IONIZED = 2
ABL_IONIZ_MODE_SELF_CONSISTENT = 3

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

MAGNETIC_FIELD_DEPENDENCE_MODE_NEGLECT = 1
MAGNETIC_FIELD_DEPENDENCE_MODE_JOREK = 2

SHIFT_MODE_NEGLECT=1
SHIFT_MODE_PRESCRIBED=2
SHIFT_MODE_ANALYTICAL=3

ZMolarMassList=[1,1,10]
isotopesMolarMassList=[2,0,0]# 0 means naturally occuring mix
molarMassList=[0.0020141,0.001008,0.020183]# kg/mol

ZSolidDensityList=[1,1,10]
isotopesSolidDensityList=[2,0,0]
solidDensityList=[205.9,86,1444]# kg/m^3

class SPI(UnknownQuantity):
    

    def __init__(self, settings, rp=None, vp=None, xp=None, t_delay = None, VpVolNormFactor=1, rclPrescribedConstant=0.01, velocity=VELOCITY_MODE_NONE, ablation=ABLATION_MODE_NEGLECT, deposition=DEPOSITION_MODE_NEGLECT, heatAbsorbtion=HEAT_ABSORBTION_MODE_NEGLECT, cloudRadiusMode=CLOUD_RADIUS_MODE_NEGLECT, magneticFieldDependenceMode=MAGNETIC_FIELD_DEPENDENCE_MODE_NEGLECT, abl_ioniz=ABL_IONIZ_MODE_NEUTRAL, shiftMode = 
SHIFT_MODE_NEGLECT, TDrift = None, T0Drift = None, DeltaYDrift = None, RmDrift = None, ZavgDriftArray = None, ZsDrift = None, isotopesDrift = None):
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
        :param int magneticFieldDependenceMode: Mode used for calculating the magnetic field dependence of the albation
        :param int shift: Model used for determining the cloud drift
        """
        super().__init__(settings=settings)

        self.velocity                    = int(velocity)
        self.ablation                    = int(ablation)
        self.deposition                  = int(deposition)
        self.heatAbsorbtion              = int(heatAbsorbtion)
        self.cloudRadiusMode             = int(cloudRadiusMode)
        self.VpVolNormFactor             = VpVolNormFactor
        self.rclPrescribedConstant       = rclPrescribedConstant
        self.magneticFieldDependenceMode = int(magneticFieldDependenceMode)
        self.abl_ioniz                   = int(abl_ioniz)
        self.shift                       = int(shiftMode)

        self.rp       = None
        self.vp       = None
        self.xp       = None
        self.t_delay  = None 
        self.nbrShiftGridCell = None

        self.TDrift        = None
        self.T0Drift       = 0
        self.DeltaYDrift  = 0
        self.RmDrift       = -1
        self.ZavgDriftArray= [0.]
        self.ZsDrift       = [0]
        self.isotopesDrift = [0]


    def setInitialData(self, rp=None, vp=None, xp=None, t_delay=None, nbrShiftGridCell = None, TDrift = None):

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
            
        if t_delay is not None:
            if np.isscalar(t_delay):
                self.t_delay = np.asarray([t_delay])
            else: self.t_delay = np.asarray(t_delay)
            
        if nbrShiftGridCell is not None:
            if np.isscalar(nbrShiftGridCell):
                self.nbrShiftGridCell = np.asarray([nbrShiftGridCell])
            else: self.nbrShiftGridCell = np.asarray(nbrShiftGridCell)

        if TDrift is not None:
            if np.isscalar(TDrift):
                self.TDrift = np.asarray([TDrift])
            else: self.TDrift = np.asarray(TDrift)

    def rpDistrParksStatistical(self,rp,kp):
        """
        Evaluates the shard size distribution function referred to as the 
        'statistical model' in P. Parks 2016 GA report (DOI:10.2172/1344852)
        """
        return kn(0,rp*kp)*kp**2*rp
        
    def sampleRpDistrParksStatistical(self, N, kp, random=np.random):
        """
        Samples N shard radii according to the distribution function 
        given by rpDistrParksStatistical()
        """
        # First we calculate the cdf, and then interpolate the cdf-values 
        # back to the corresponding radii at N randomly chosen points between 0 and 1
        rp_integrate=np.linspace(1e-10/kp,10/kp,5000)
        cdf=integrate.cumulative_trapezoid(y=self.rpDistrParksStatistical(rp_integrate,kp),x=rp_integrate)
        return np.interp(random.uniform(size=N),np.hstack((0.0,cdf)),rp_integrate)
        
    def setRpParksStatistical(
        self, nShard, Ninj, Zs, isotopes, molarFractions, ionNames,
        opacity_modes=None, add=True, n=1e0, random=np.random,
        charged_advection_modes=None, charged_prescribed_advections=None,
        rChargedPrescribedAdvections=None, tChargedPrescribedAdvections=None,
        neutral_advection_modes=None, neutral_prescribed_advections=None,
        rNeutralPrescribedAdvections=None, tNeutralPrescribedAdvections=None,
        charged_diffusion_modes=None, charged_prescribed_diffusions=None,
        rChargedPrescribedDiffusions=None, tChargedPrescribedDiffusions=None,
        neutral_diffusion_modes=None, neutral_prescribed_diffusions=None,
        rNeutralPrescribedDiffusions=None, tNeutralPrescribedDiffusions=None,
        **kwargs
    ):
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
        counter=0
        for iZ in range(len(Zs)):
            for iList in range(len(solidDensityList)):
                if Zs[iZ]==ZSolidDensityList[iList] and isotopes[iZ]==isotopesSolidDensityList[iList]:
                    solidDensityIZ=solidDensityList[iList]
                if Zs[iZ]==ZMolarMassList[iList] and isotopes[iZ]==isotopesMolarMassList[iList]:
                    molarMassIZ=molarMassList[iList]
                else:
                    counter+=1
            if counter==len(solidDensityList):
                raise EquationException("spi: Pellet type is not recognized. Currently only neon and deuterium pellets are supported. To support other types fill in the material data in src/Equations/SPIHandler.cpp and py/DREAM/Settings/Equations/SPI.py")
            counter=0
            molarVolume+=molarFractions[iZ]*molarMassIZ/solidDensityIZ
            
        solidParticleDensity=N_A/molarVolume
       
       
        # Calculate inverse characteristic shard size
        kp=(6*np.pi**2*solidParticleDensity*nShard/Ninj)**(1/3)
        
        # Sample the shard sizes and rescale to get exactly the 
        # specified number of particles in the pellet
        rp_init=self.sampleRpDistrParksStatistical(nShard, kp, random=random)
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
            
        # Fix arrays of settings to have correct shape if None is specified        
        if opacity_modes is None:
            opacity_modes = []
            for iZ in range(len(Zs)):
                opacity_modes.append(Ions.ION_OPACITY_MODE_TRANSPARENT)
                
        if charged_advection_modes is None:
            charged_advection_modes = []
            for iZ in range(len(Zs)):
                charged_advection_modes.append(Ions.ION_CHARGED_ADVECTION_MODE_NONE)
        if charged_prescribed_advections is None:
            charged_prescribed_advections = []
            for iZ in range(len(Zs)):
                charged_prescribed_advections.append(None)
        if rChargedPrescribedAdvections is None:
            rChargedPrescribedAdvections = []
            for iZ in range(len(Zs)):
                rChargedPrescribedAdvections.append(None)
        if tChargedPrescribedAdvections is None:
            tChargedPrescribedAdvections = []
            for iZ in range(len(Zs)):
                tChargedPrescribedAdvections.append(None)
                
        if neutral_advection_modes is None:
            neutral_advection_modes = []
            for iZ in range(len(Zs)):
                neutral_advection_modes.append(Ions.ION_NEUTRAL_ADVECTION_MODE_NONE)
        if neutral_prescribed_advections is None:
            neutral_prescribed_advections = []
            for iZ in range(len(Zs)):
                neutral_prescribed_advections.append(None)
        if rNeutralPrescribedAdvections is None:
            rNeutralPrescribedAdvections = []
            for iZ in range(len(Zs)):
                rNeutralPrescribedAdvections.append(None)
        if tNeutralPrescribedAdvections is None:
            tNeutralPrescribedAdvections = []
            for iZ in range(len(Zs)):
                tNeutralPrescribedAdvections.append(None)
                
        if charged_diffusion_modes is None:
            charged_diffusion_modes = []
            for iZ in range(len(Zs)):
                charged_diffusion_modes.append(Ions.ION_CHARGED_DIFFUSION_MODE_NONE)
        if charged_prescribed_diffusions is None:
            charged_prescribed_diffusions = []
            for iZ in range(len(Zs)):
                charged_prescribed_diffusions.append(None)
        if rChargedPrescribedDiffusions is None:
            rChargedPrescribedDiffusions = []
            for iZ in range(len(Zs)):
                rChargedPrescribedDiffusions.append(None)
        if tChargedPrescribedDiffusions is None:
            tChargedPrescribedDiffusions = []
            for iZ in range(len(Zs)):
                tChargedPrescribedDiffusions.append(None)
                
        if neutral_diffusion_modes is None:
            neutral_diffusion_modes = []
            for iZ in range(len(Zs)):
                neutral_diffusion_modes.append(Ions.ION_NEUTRAL_DIFFUSION_MODE_NONE)
        if neutral_prescribed_diffusions is None:
            neutral_prescribed_diffusions = []
            for iZ in range(len(Zs)):
                neutral_prescribed_diffusions.append(None)
        if rNeutralPrescribedDiffusions is None:
            rNeutralPrescribedDiffusions = []
            for iZ in range(len(Zs)):
                rNeutralPrescribedDiffusions.append(None)
        if tNeutralPrescribedDiffusions is None:
            tNeutralPrescribedDiffusions = []
            for iZ in range(len(Zs)):
                tNeutralPrescribedDiffusions.append(None)
                                
        # Add an ion species connected to this pellet to the ion settings
        for iZ in range(len(Zs)):
            
            # SPIMolarFraction must have the smae length as all pellet shard, 
            # not only the pellet which is initiated here, so set the molar fraction 
            # to zero for previously set shards
            SPIMolarFraction=np.zeros(len(self.rp))
            SPIMolarFraction[-nShard:]=molarFractions[iZ]*np.ones(nShard)
            
            self.settings.eqsys.n_i.addIon(
                name=ionNames[iZ], n=n, Z=Zs[iZ], isotope=isotopes[iZ], opacity_mode=opacity_modes[iZ], iontype=Ions.IONS_DYNAMIC_NEUTRAL,
                SPIMolarFraction=SPIMolarFraction, charged_diffusion_mode = charged_diffusion_modes[iZ],
                charged_prescribed_diffusion = charged_prescribed_diffusions[iZ], rChargedPrescribedDiffusion = rChargedPrescribedDiffusions[iZ],
                tChargedPrescribedDiffusion = tChargedPrescribedDiffusions[iZ], neutral_diffusion_mode = neutral_diffusion_modes[iZ],
                neutral_prescribed_diffusion = neutral_prescribed_diffusions[iZ], rNeutralPrescribedDiffusion = rNeutralPrescribedDiffusions[iZ],
                tNeutralPrescribedDiffusion = tNeutralPrescribedDiffusions[iZ], charged_advection_mode = charged_advection_modes[iZ],
                charged_prescribed_advection = charged_prescribed_advections[iZ], rChargedPrescribedAdvection = rChargedPrescribedAdvections[iZ],
                tChargedPrescribedAdvection = tChargedPrescribedAdvections[iZ], neutral_advection_mode = neutral_advection_modes[iZ],
                neutral_prescribed_advection = neutral_prescribed_advections[iZ], rNeutralPrescribedAdvection = rNeutralPrescribedAdvections[iZ],
                tNeutralPrescribedAdvection = tNeutralPrescribedAdvections[iZ], **kwargs
            )
            
              
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
            
    def setShardVelocitiesUniform(
        self, nShard, abs_vp_mean, abs_vp_diff, alpha_max, tilt=0,
        t_delay = 0, nDim=2,add=True, shards=None, random=np.random
    ):
        """
        Sets self.vp to a vector storing the (x,y,z)-components of nShard shard velosities,
        assuming a uniform velocity distribution over a nDim-dimensional cone whose axis
        is anti-parallell to the x-axis. TODO: implement support for an arbitrary axis?
        
        :param int nShard: Number of shards
        :param float abs_vp_mean: Mean of the magnitude of the shard velocities
        :param float abs_vp_diff: width of the uniform distribution of the magnitude of the shard velocities
        :param float alpha_max: Span of divergence angle (ie twice the opening angle of the cone)
        :param float tilt: Tilt angle of the SPI w.r.t. the x-axis (in radians; positive angle is up).
        :param int nDim: number of dimensions into which the shards should be spread
        :param bool add: If 'True', add the new pellet shard velocities to the existing ones, otherwise 
             existing shards are cleared
        :param slice shards: indices of existing shards whose velocities should be updated. If not 'None', 
                add is set to 'False' and nShard is set to the number of indices to be updated
        :param random: Random number generator to use (default: numpy.random).
        """
        
        if shards is not None:
        	nShard=len(self.vp[shards])
        	add=False
        	
        if np.isscalar(t_delay): 
            t_delay = t_delay*np.ones(nShard)
        
        # Sample magnitude of velocities
        abs_vp_init=(abs_vp_mean+abs_vp_diff*(-1+2*random.uniform(size=nShard)))
        
        # Sample directions uniformly over a nDim-dimensional cone and set the velocity vectors
        vp_init=np.zeros(3*nShard)
        if nDim==1:
            # in 1D, the "cone" simply becomes a straight line
            vp_init[0::3]=-abs_vp_init
            
        elif nDim==2:
            # in 2D, the cone becomes a circle sector
            alpha=alpha_max*(-1+2*random.uniform(size=nShard)) + tilt
            vp_init[0::3]=-abs_vp_init*np.cos(alpha)
            vp_init[1::3]=abs_vp_init*np.sin(alpha)
            
        elif nDim==3:
            # The solid angle covered by the part of the cone between alpa and d(alpha) 
            # is proportional to sin(alpha), and the normalised probability distribution 
            # becomes f(alpha)=sin(alpha)/(1-cos(alpha_max/2)). We sample from this
            # distribution by applying the inverse cdf to uniformly drawn numbers
            # between 0 and 1
            alpha=np.arccos(1-random.uniform(size=nShard)*(1-np.cos(alpha_max/2)))
            
            # The angle in the yz-plane is simply drawn randomly
            phi=2*np.pi*random.uniform(size=nShard)
            
            # Finally calculate the velocity vectors
            vp_init[0::3]=-abs_vp_init*np.cos(alpha)
            vp_init[1::3]=abs_vp_init*np.sin(alpha)*np.cos(phi)
            vp_init[2::3]=abs_vp_init*np.sin(alpha)*np.sin(phi)
            
        else:
            raise EquationException("spi: Invalid number of dimensions into which the pellet shards are spread")
            
        if add and self.vp is not None:
            self.vp=np.concatenate((self.vp,vp_init))
            self.t_delay=np.concatenate((self.t_delay,t_delay))
        elif shards is not None:
        	# Pick out the components of the stored shard velocities...
        	vpx=self.vp[0::3]
        	vpy=self.vp[1::3]
        	vpz=self.vp[2::3]
        	
        	# ... Change the velocities of the shards specified in the input...
        	vpx[shards]=vp_init[0::3]
        	vpy[shards]=vp_init[1::3]
        	vpz[shards]=vp_init[2::3]
        	
        	# ...and finally set the stored velocities to the updated ones
        	self.vp[0::3]=vpx
        	self.vp[1::3]=vpy
        	self.vp[2::3]=vpz
        	
        	self.t_delay[shards] = t_delay
        else:
            self.vp=vp_init
            self.t_delay = t_delay
            
    def setParamsVallhagenMSc(
        self, nShard, Ninj, Zs, isotopes, molarFractions, ionNames,
        shatterPoint, abs_vp_mean,abs_vp_diff,alpha_max,t_delay=0,
        tilt=0, nDim=2, add=True, opacity_modes=None, nbrShiftGridCell=0,
        TDrift=None, random=np.random, **kwargs
    ):
        """
        Wrapper for setRpParksStatistical(), setShardPositionSinglePoint() and setShardVelocitiesUniform(),
        which combined are used to set up an SPI-scenario similar to those in Oskar Vallhagens MSc thesis
        (available at https://hdl.handle.net/20.500.12380/302296).
        """
        
        kp=self.setRpParksStatistical(nShard, Ninj, Zs, isotopes, molarFractions, ionNames, opacity_modes, add, random=random, **kwargs)
        self.setShardPositionSinglePoint(nShard,shatterPoint,add)
        self.setShardVelocitiesUniform(
            nShard=nShard, abs_vp_mean=abs_vp_mean, abs_vp_diff=abs_vp_diff,
            alpha_max=alpha_max, tilt=tilt, t_delay=t_delay, nDim=nDim, add=add,
            random=random
        )
        
        if add and self.nbrShiftGridCell is not None:
            self.nbrShiftGridCell = np.concatenate((self.nbrShiftGridCell,nbrShiftGridCell*np.ones(nShard, dtype=np.int64)))
        else:
            self.nbrShiftGridCell = nbrShiftGridCell*np.ones(nShard, dtype=np.int64)
            
        # Perhaps it would be better to force the user to explicitly set the shift mode,
        # but this helps to ensure backwards compatibility with scripts relying on that 
        # setting nbrShiftGridCell>0 automatically gives a corresponding prescribed shift
        if nbrShiftGridCell>0:
            self.setShift(SHIFT_MODE_PRESCRIBED)
            
        if TDrift is not None:
            if np.isscalar(TDrift):
                TDrift = TDrift*np.ones(nShard)
            if add and self.TDrift is not None:
                self.TDrift = np.concatenate((self.TDrift,TDrift))
            else:
                self.TDrift = TDrift

        return kp
        
    def setShiftParamsPrescribed(self, shift = SHIFT_MODE_PRESCRIBED, nbrShiftGridCell=None, add=True):
        self.setShift(shift)
        if nbrShiftGridCell is not None:
            if add and self.nbrShiftGridCell is not None:
                self.nbrShiftGridCell = np.concatenate((self.nbrShiftGridCell,nbrShiftGridCell))
            else:
                self.nbrShiftGridCell = nbrShiftGridCell
        
    def setShiftParamsAnalytical(self, shift = SHIFT_MODE_ANALYTICAL, TDrift=None, T0Drift=0, DeltaYDrift=0, RmDrift=-1, ZavgDriftArray=[0.], ZsDrift=[0], isotopesDrift=[0], add=True):
        """
        Specifies model parameters to be used for calculating the shift. Apart from the shift mode-argument, the parameters below apply to SHIFT_MODE_ANALYTICAL
        
        :param int shift: Model used for determining the cloud drift
        :param float T0Drift: cloud temperature close to the pellet (before the cloud has drifted away from the pellet)
        :param numpy.ndarray TDrift: representative cloud temperature during the drift motion for each shard
        :param float DeltaYDrift: characteristic half-thickness of the drifting cloud (should be similar to the radius of the neutral cloud around the pellet)
        :param float RmDrift: major radius of the magnetic axis, only used if the major radius is otherwise infinite in the simulation
        :param list ZavgDriftArray: average charge states inside the drifting cloud of all drifting ion species. These can not be calculated using the ADAS rates because the conditions in the drifting cloud, especially the density and optical thickness, are very different from the validity range and assumptions in ADAS, and we therefore take user-given estimates for them. Note that his list does NOT neccessarily have the same shape as the list of atomic numbers and isotopes included in the simulation, but instead the ZavgDriftArray-list and the ZsDrift and isotopesDrift-lists below will instead be used to look up the average charge state inside the drifting cloud for all the simulated ion species included in the pellet.
        :param list ZsDrift: atomic numbers of all the drifting ion species, corresponding to the average charge states listed in the ZavgDriftArray-list above
        :param list isotopesDrift: isotopes of all the drifting ion species, corresponding to the average charge states listed in the ZavgDriftArray-list above
        """
        self.setShift(shift)
        if TDrift is not None:
            if add and self.TDrift is not None:
                if np.isscalar(TDrift):
                    self.TDrift = np.concatenate((self.TDrift, [TDrift]))
                else:
                    self.TDrift = np.concatenate((self.TDrift, TDrift))
            else:
                self.TDrift = TDrift
        self.T0Drift = T0Drift
        self.DeltaYDrift = DeltaYDrift
        self.RmDrift = RmDrift
        self.ZavgDriftArray = ZavgDriftArray
        self.ZsDrift = ZsDrift
        self.isotopesDrift = isotopesDrift
        
    def setVpVolNormFactor(self,VpVolNormFactor):
        self.VpVolNormFactor=VpVolNormFactor

    def setRclPrescribedConstant(self,rclPrescribedConstant):
        self.rclPrescribedConstant=rclPrescribedConstant
        
    def shiftTimeDelay(self, tShift):
        self.t_delay = self.t_delay - tShift


    def resetTimeDelay(self):
        self.t_delay = 0


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
        
    def setMagneticFieldDependenceMode(self, magneticFieldDependenceMode):
        """
        Specifies which model to use for calculating the
        magnetic field dependence of the ablation
        """
        self.magneticFieldDependenceMode = int(magneticFieldDependenceMode)
        
    def setAblIoniz(self, abl_ioniz):
        """
        Specifies which model to use for calculating the
        charge state distribution with which the recently 
        ablated material is deposited
        """
        self.abl_ioniz = int(abl_ioniz)
        
    def setShift(self, shift):
        """
        Specifies which model to use for calculating the
        ablation cloud drift
        """
        self.shift = int(shift)

    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        if 'velocity' in data:
            self.velocity       = int(data['velocity'])
        if 'ablation' in data:
            self.ablation       = int(data['ablation'])
        if 'deposition' in data:
            self.deposition     = int(data['deposition'])
        if 'shift' in data:
            self.shift          = int(data['shift'])
        if 'TDrift' in data:
            self.TDrift              = [float(x) for x in data['TDrift']]
        if 'T0Drift' in data:
            self.T0Drift             = float(data['T0Drift'])
        if 'DeltaYDrift' in data:
            self.DeltaYDrift        = float(data['DeltaYDrift'])
        if 'RmDrift' in data:
            self.RmDrift             = float(data['RmDrift'])
        if 'ZavgDriftArray' in data:
            self.ZavgDriftArray      = [float(x) for x in data['ZavgDriftArray']]
        if 'ZsDrift' in data:        
            self.ZsDrift             = [float(x) for x in data['ZsDrift']]
        if 'isotopesDrift' in data:        
            self.isotopesDrift             = [float(x) for x in data['isotopesDrift']]
        if 'heatAbsorbtion' in data:
            self.heatAbsorbtion = int(data['heatAbsorbtion'])
        if 'cloudRadiusMode' in data:
            self.cloudRadiusMode = int(data['cloudRadiusMode'])
        if 'magneticFieldDependenceMode' in data:
            self.magneticFieldDependenceMode = int(data['magneticFieldDependenceMode'])
        if 'abl_ioniz' in data:
            self.abl_ioniz = int(data['abl_ioniz'])
            

        if 'VpVolNormFactor' in data:
            self.VpVolNormFactor = data['VpVolNormFactor']
        if 'rclPrescribedConstant' in data:
            self.rclPrescribedConstant = data['rclPrescribedConstant']
        if 'nbrShiftGridCell' in data:
            self.nbrShiftGridCell = data['nbrShiftGridCell']

        if 'init' in data:
            if 'rp' in data['init']:
                self.rp              = data['init']['rp']
            if 'vp' in data['init']:
                self.vp              = data['init']['vp']
            if 'xp' in data['init']:
                self.xp              = data['init']['xp']
            if 't_delay' in data['init']:
                self.t_delay         = data['init']['t_delay']


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this SPI object.
        """
        # If no SPI settings have been given, set everything to zero (to avoid a DREAMIOException)
        # Before this stage it is usefull to use None to indicate if any SPI settings have been made yet,
        # to know if there are any previous shards to add the new ones to, so therefore
        # we don't set this default setting until this stage
        if self.t_delay is None:
            if self.rp is not None:
                self.t_delay=np.zeros(self.rp.shape)
            else:
                self.t_delay=np.array([0])

        data = {
            'velocity': self.velocity,
            'ablation': self.ablation,
            'deposition': self.deposition,
            'shift': self.shift,
            'T0Drift': self.T0Drift,
            'DeltaYDrift': self.DeltaYDrift,
            'RmDrift': self.RmDrift,
            'ZavgDriftArray': self.ZavgDriftArray,
            'ZsDrift': self.ZsDrift,
            'isotopesDrift': self.isotopesDrift,
            'heatAbsorbtion': self.heatAbsorbtion,
            'cloudRadiusMode': self.cloudRadiusMode,
            'magneticFieldDependenceMode': self.magneticFieldDependenceMode,
            'abl_ioniz': self.abl_ioniz,
            'VpVolNormFactor': self.VpVolNormFactor,
            'rclPrescribedConstant': self.rclPrescribedConstant
        }
        
            
        data['init'] = {}
        
        if self.rp is not None:
            data['init']['rp']=self.rp
        if self.vp is not None:
            data['init']['vp']=self.vp
        if self.xp is not None:    
            data['init']['xp']=self.xp
        if self.t_delay is not None: 
            data['init']['t_delay']=self.t_delay
            
        if self.nbrShiftGridCell is not None:
            data['nbrShiftGridCell'] = self.nbrShiftGridCell
        """
        else:
            if self.rp is not None:
                self.nbrShiftGridCell = np.zeros(self.rp.shape)
            else:
                self.nbrShiftGridCell = np.array([0])
        """

        if self.TDrift is not None:
            data['TDrift'] = self.TDrift
        """
        else:
            if self.rp is not None:
                self.TDrift = np.zeros(self.rp.shape)
            else:
                self.TDrift=np.array([0])
        """

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
        if type(self.shift) != int:
            raise EquationException("spi: Invalid value assigned to 'shift'. Expected integer.")
        if self.shift == SHIFT_MODE_ANALYTICAL:
            if self.T0Drift<0: 
                raise EquationException("spi: Invalid value assigned to 'T0Drift'. Expected positive float.")
            if any(self.TDrift)<=0:
                raise EquationException("spi: Invalid value assigned to 'TDrift'. Expected array of positive floats.")
            if self.DeltaYDrift<0:
                raise EquationException("spi: Invalid value assigned to 'DeltaYDrift'. Expected positive float.")
            if self.RmDrift<0 and self.RmDrift!=-1:
                raise EquationException("spi: Invalid value assigned to 'RmDrift'. Expected positive float.")
            if len(self.ZavgDriftArray)!=len(self.ZsDrift):
                raise EquationException("spi: Invalid value assigned to 'ZavgDriftArray'. Expected array of positive floats with the same shape as 'ZsDrift'.")
            if len(self.isotopesDrift)!=len(self.ZsDrift):
                raise EquationException("spi: Invalid value assigned to 'isotopesDrift'. Expected array of positive floats with the same shape as 'ZsDrift'.")
            if self.deposition!=DEPOSITION_MODE_LOCAL:
                raise EquationException("spi: Invalid value assigned to 'shift'. To enable shift activate deposition.")
        if type(self.heatAbsorbtion) != int:
            raise EquationException("spi: Invalid value assigned to 'heatAbsorbtion'. Expected integer.")



    def verifySettingsPrescribedInitialData(self):
        if vp.size!=3*rp.size:
            raise EquationException("Missmatch in size of initial data arrays for rp and vp. Expected vp to have a size 3 times the size of rp")
        if xp.size!=3*rp.size:
            raise EquationException("Missmatch in size of initial data arrays for rp and xp. Expected xp to have a size 3 times the size of rp")
        
