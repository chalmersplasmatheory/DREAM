#!/usr/bin/env python3
#
# This script shows how to set up an D/Ne SPI scenario in an ITER-like setting,
# using the ability to read initial values from CORSICA data with the ITER_profile_reader-tool,
# numerical geometry in the LUKE-format (converted from GQDESK-data), adn electron transport
# coefficients with the transport_coeffs_reader-tool, depending on the settings.
# The injection can either be separated into one stage with pure D and one stage with pure NE,
# or be made as a single stage injection with a similar total amount of particles.
#
################################################################################################ 

import numpy as np
import scipy as sp
import sys
import h5py

from scipy import integrate
from scipy.special import kn

sys.path.append('../../py/')
sys.path.append('../../tools/')

from DREAM.DREAMSettings import DREAMSettings
from DREAM.DREAMOutput import DREAMOutput
from DREAM import runiface
from ITER_Profile_reader import ITER_Profile_reader
from transport_coeffs_reader import transport_coeffs_reader
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.Ions as IonsAll
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ElectricField as Efield
import DREAM.Settings.Equations.RunawayElectrons as RE
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.TimeStepper as TimeStep
import DREAM.Settings.Equations.SPI as SPI
import DREAM.Settings.TransportSettings as Transport
import DREAM.Settings.RadialGrid as RGrid


from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature

# Makes sure this script generates exactly the same result every time the same settings are used
np.random.seed(1)

ds = DREAMSettings()

#####################
# Scenario settings #
#####################

# Specify scenario, folder with scenario data and what time stamp to take the data from
scenario_name = 'H26'
scenario_folder = 'H26'
scenario_time = '00060'

#scenario_name = 'DTHMODE24'
#scenario_folder = 'DThmode24'
#scenario_time = '00400'

# Initial position (shatter point) for all the pellet shards (R, Z)
shatterPoint = np.array([8.568, 0.6855])

# Settings for the first SPI (presumably consisting mostly of deuterium)
nShardD=478 # Number of shards
NinjD=1.85e24 # Number of atoms
alpha_maxD=0.17 # Divergence angle
abs_vp_meanD=500 # Mean shard speed
abs_vp_diffD=0.4*abs_vp_meanD # Width of the uniform shard speed distribution
molarFractionNe=5/185 # Molar fraction of neon (the rest is deuterium)
# molarFractionNe = 0

# Settings for the second Neon SPI
nShardNe=0
NinjNe=0
alpha_maxNe=0.17
abs_vp_meanNe=200
abs_vp_diffNe=0.4*abs_vp_meanNe

# Specify data file to read electron transport coefficients from, if any
#transp_coeffs_filename = 'transp_coeffs_Di_7_5_MA'
transp_coeffs_filename = None

# If no file with transport coefficients is specified, 
# set homogenious transport coefficients to give the desired TQ time t_duration
# The appropriate ratio of t_duration_over_t_diffusion is determined from parameter scans
# presented at the ITER meeting from 2022-07-22
if transp_coeffs_filename is None:
    t_duration = 3e-3 # s
    if scenario_name == 'H26':
        t_duration_over_t_diffusion = 35
    elif scenario_name == 'DTHMODE24':
        t_duration_over_t_diffusion = 65
else:
    t_duration = None
    t_duration_over_t_diffusion = None
    
# Set parameters for exponentially decaying ion transport coeffcicients
t_exp_ion = 0.5e-3
DrrIon = 4e3
ArIon = -2e3
  
# TQ onset criterion, see presentation from the ITER meeting 2022-07-22 for further description  
TQ_onset = 'Late'
# TQ_onset = 'Early'
# TQ_onset = None # Use this to run a simulation without transport that can be used to determine when the TQ should start

if TQ_onset is None:
    use_heat_transport=False
    use_nre_transport=False
    use_ion_transport = False
else:
    use_heat_transport=True
    use_nre_transport=True
    use_ion_transport = True
    
# If assuming a late TQ, specify the file (filename_output) of the simulation 
# (ran without transport coefficients) to be used to determine the TQ onset time
if TQ_onset == 'Late':
    folder_output = 'output/'
    filename_output='output_restart_injection_ITER_'+scenario_name+'_INDEX_SPI_nShardD'+str(nShardD)+'NNe'+str(NinjD*molarFractionNe)+'abl_ioniz1radius_wall3new.h5'
    
#######################################
# Physical model parameters  #
#######################################

# Runaway settings
hotTailGrid_enabled = False
use_fluid_runaways = True

# Set fluid runaway generation rates
if use_fluid_runaways:
	if scenario_name == 'DTHMODE24':
		ds.eqsys.n_re.setCompton(RE.COMPTON_RATE_ITER_DMS)
		ds.eqsys.n_re.setTritium(True)
	ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID_HESSLOW)
	ds.eqsys.n_re.setDreicer(RE.DREICER_RATE_NEURAL_NETWORK)
	if not hotTailGrid_enabled:
		ds.eqsys.n_re.setHottail(RE.HOTTAIL_MODE_ANALYTIC_ALT_PC)
	
# Settings for the Svensson fluid runaway transport	
ds.eqsys.n_re.transport.type = Transport.TRANSPORT_SVENSSON
pstar = 0.3
ds.eqsys.n_re.transport.setSvenssonPstar(pstar)

# Kinetic grid settings
if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# set collision settings
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
#ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_NEGLECT
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
#ds.collisions.lnlambda = Collisions.LNLAMBDA_CONSTANT
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL



# Set models for advancing the shard motion and ablation
ds.eqsys.spi.setVelocity(SPI.VELOCITY_MODE_PRESCRIBED) # Constant prescribed velocities
ds.eqsys.spi.setAblation(SPI.ABLATION_MODE_FLUID_NGS) # Parks NGS formula based on T_cold
# ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL_GAUSSIAN) # Use a gaussian deposition kernel of width rcl
ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL) # Delta function deposition kernel
# ds.eqsys.spi.setHeatAbsorbtion(SPI.HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS) # Remove all heat flowing through a disc 
#                                                                            of radius rcl from background plasma.
#                                                                            when assuming a local and immediate deposition
#                                                                            the heat absorbed by the ablated material is
#                                                                            immediately returned, so this term should then 
#                                                                            be switched off 
# ds.eqsys.spi.setCloudRadiusMode(SPI.CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT)
ds.eqsys.spi.setMagneticFieldDependenceMode(SPI.MAGNETIC_FIELD_DEPENDENCE_MODE_JOREK)
ds.eqsys.spi.setAblIoniz(SPI.ABL_IONIZ_MODE_NEUTRAL)

# Size of the neutral cloud (only relevant when using gaussian deposition or heat absorbtion)
rcl=0.01
ds.eqsys.spi.setRclPrescribedConstant(rcl)

#######################################
# Numerical simulation parameters #
#######################################

# Choose which parts of the disruption should be simulated
run_init=True # Includes a single-step run to extract the conductivity and 
#                another one to set the efield according to the wanted current profiel

run_injection_init=True # Accelerate the distribution function to carry the right current

run_injection=True # D or D/Ne injection
run_CQ=True # Second Ne injection (if any), beginning of the CQ

# Specify number of restarts to do during the CQ
nCQ_restart_start=2 # Number of CQ restart to start from (usefull for re-running only a later part of the CQ)
nCQ_restart=1 # How many CQ restarts to run

# Time steps during the various restarts
Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps

Tmax_init2 = 1e-3   
Nt_init2 = 300      

# For first stage of a staggered injection
#Tmax_injection = 3.4e-3
#Nt_injection = 1360   

# For single stage
Tmax_injection = 6e-3
Nt_injection = 12000  

# For the start of the CQ (second injection for a staggered injection)
Tmax_CQ = 17e-3
Nt_CQ = 30000

# For the first (currently only) restart during the CQ
Tmax_CQ_restart = 179.6e-3
Nt_CQ_restart = 40000

# Grid resolution parameters
Nr = 11             # number of radial grid points
Np = 120            # number of momentum grid points
Nxi = 5             # number of pitch grid points
pMax = 3            # maximum momentum in m_e*c

times  = [0]        # times at which parameters are given

# Set up radial grid
ds.radialgrid.setType(RGrid.TYPE_NUMERICAL)
ds.radialgrid.setNumerical('ITER_'+scenario_name+'_LUKE.h5', format=RGrid.FILE_FORMAT_LUKE)
# radius_wall = 2.15  # location of the wall 
radius_wall = 3  # location of the wall 
ds.radialgrid.setWallRadius(radius_wall)
ds.radialgrid.setNr(Nr)

# Extract further data for the grid and set accordingly (do you still need to set major and minor radius manually like this?)
dataLUKE = h5py.File('ITER_'+scenario_name+'_LUKE.h5','r')['equil']
ds.radialgrid.setMajorRadius(np.array([dataLUKE['Rp']]))
ds.radialgrid.setMinorRadius(ds.radialgrid.a-1e-6) # To avoid interpolation errors at the edge

# Set advection interpolation method for the ion transport
ds.eqsys.n_i.setAdvectionInterpolationMethodCharged(ad_int=IonsAll.AD_INTERP_UPWIND,
        ad_jac=IonsAll.AD_INTERP_JACOBIAN_UPWIND, fluxlimiterdamping=1.0)
  
# Boundary conditions for the thermal energy transport    
ds.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_F_0)

# Boundary condition for the Svensson fluid runaway transport
ds.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
        
# Hot-tail grid settings (only relevant if the hot-tail grid is enabled)
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)
ds.hottailgrid.setBiuniformGrid(psep=0.07,npsep=70)
ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
#ds.eqsys.f_hot.enableIonJacobian(False)

# Settings for the nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_MKL)
ds.solver.setMaxIterations(maxiter = 500)
ds.solver.setTolerance(reltol=0.001)
#ds.solver.setVerbose(True)

#################################
# Set initial plasma parameters #
#################################

# The current is set later after the first initalisation run, 
# as one needs flux surface averages calculated during runtime 
# to rescale the current data from the scenario data file 
# to the current density at the outboard midplane which should be given as input to DREAM

# Create an ITER_profile_reader-object containing the initial data from the scenario data file
I = ITER_Profile_reader(scenario_folder+'/P_'+scenario_name+'_ITER_MR_'+scenario_time+'.TXT')

# Minor radius coordinates at which the input profiles are given
rref = I.load_profile('R_minor')

# Set E_field 
E_initial = 0.00032   # initial electric field in V/m (will be overwritten later when the current is set)
E_wall = 0.0          # boundary electric field in V/m
inverse_wall_time = 2 # resistive time scale of the wall in s^{-1}
I0_add_spike = 0.0 # Current to add to the initial current from the scenario data file to emulate the current spike
# set prescribed and flat E-field profile and defualt boundary conditions for now (will be changed after the initialisation runs)
efield = E_initial*np.ones((len(times), len(rref))) 
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=rref)
ds.eqsys.E_field.setBoundaryCondition()

# Set temperature profile
temperature = I.load_profile('Te')
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=rref)

# Set background ions
n_D = I.load_profile('Ni')
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=0.5*n_D, r=rref, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
ds.eqsys.n_i.addIon(name='T', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=0.5*n_D, r=rref, tritium=True, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)

# Initial distribution function (needed also for fluid hot-tail)
nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=rref, T0=temperature.flatten())

#######################################
# Set electron transport coefficients #
#######################################

# Ion transport coefficients will be set later as the SPI settings 
# and the settings for the ion species corresponding to the SPI is most convenient to set simultaneously

if use_nre_transport or use_heat_transport:
    # File to find TQ onset criterion from for late TQ, if any
    if TQ_onset == 'Late':
        t_transp_onset = transport_coeffs_reader.tBeforeOnsetFromQCritAndTcoldFromOutput(q=I.load_profile('q'), qcrit=2, rref=rref, Tcrit = 10, filename=folder_output+filename_output)
    elif TQ_onset == 'Early':
        # Estimate the time it takes for the shards to reach the q=2 surface 
        # and use this as the onset time for the stochastic transport
        t_transp_onset = transport_coeffs_reader.tBeforeOnsetFromQCritAndPelletShardPosition(q=I.load_profile('q'), rref=rref, shatterPoint=shatterPoint, abs_vp_mean=abs_vp_meanD, qcrit=2)
        
    # Load transport coefficients from file if such a file is specified,
    # otherwise calculate transport coefficients to get the desired TQ time, 
    # and store in the transport_coeffs_reader object
    if transp_coeffs_filename is None:
        T = transport_coeffs_reader(transp_coeffs_filename, t_before_onset = t_transp_onset, t_duration = t_duration)
        T.setDrrFromTQTimeScaleBesselApproximation(t_duration = t_duration, a = ds.radialgrid.a, Trepr = temperature[0], t_duration_over_t_diffusion = t_duration_over_t_diffusion)
    else:
        T = transport_coeffs_reader(transp_coeffs_filename+'.h5', t_before_onset = t_transp_onset)
        
     
###########################################
# Set SPI configuration and ion transport #
###########################################

# Convert shatter point to SPI coordinates
shatterPoint = np.array([shatterPoint[0] - ds.radialgrid.R0, shatterPoint[1]-float(dataLUKE['Zp'][:]), 0])

# Input settings for the first injection stage to the settings object via the helper function 
# for set up a similar configuration as in Oskar Vallhagens MSc thesis
# The shard velocities are set to zero for now, 
# and will be changed later when the injection is supposed to start
if molarFractionNe>0:
	if use_ion_transport:
		ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardD, Ninj=NinjD, Zs=[1,10], isotopes=[2,0],
			opacity_modes=[Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE, Ions.ION_OPACITY_MODE_TRANSPARENT], 
			molarFractions=[1-molarFractionNe,molarFractionNe], ionNames=['D_inj_mix','Ne_inj_mix'], abs_vp_mean=0, abs_vp_diff=0, 
			alpha_max=alpha_maxD, shatterPoint=shatterPoint, n = np.ones(rref.size), r = rref, 
			charged_advection_modes = [Ions.ION_CHARGED_ADVECTION_MODE_PRESCRIBED, Ions.ION_CHARGED_ADVECTION_MODE_PRESCRIBED], 
	        neutral_advection_modes = [Ions.ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED, Ions.ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED], 
	        charged_diffusion_modes = [Ions.ION_CHARGED_DIFFUSION_MODE_PRESCRIBED, Ions.ION_CHARGED_DIFFUSION_MODE_PRESCRIBED], 
	        neutral_diffusion_modes = [Ions.ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED, Ions.ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED],
	        t_transp_expdecay_all_cs = t_exp_ion, t_transp_start_expdecay_all_cs = t_transp_onset, diffusion_initial_all_cs = DrrIon, advection_initial_all_cs = ArIon)
	else:
		ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardD, Ninj=NinjD, Zs=[1,10], isotopes=[2,0], opacity_modes=[Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE, Ions.ION_OPACITY_MODE_TRANSPARENT], molarFractions=[1-molarFractionNe,molarFractionNe], ionNames=['D_inj_mix','Ne_inj_mix'], abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxD, shatterPoint=shatterPoint, n = np.ones(rref.size), r = rref)
else:
	ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardD, Ninj=NinjD, Zs=[1], isotopes=[2], opacity_modes=[Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE], molarFractions=[1], ionNames=['D_inj'], abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxD, shatterPoint=shatterPoint)



# Input settings for the second injection stage (if any) to the settings object via the helper function 
# for set up a similar configuration as in Oskar Vallhagens MSc thesis
# The shard velocities are set to zero for now, 
# and will be changed later when the injection is supposed to start
if nShardNe>0:
	if use_ion_transport:
		ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardNe, Ninj=NinjNe, Zs=[10], isotopes=[0], molarFractions=[1], ionNames=['Ne_inj'], 
			abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxNe, shatterPoint=shatterPoint,   
			charged_advection_modes = [Ions.ION_CHARGED_ADVECTION_MODE_PRESCRIBED], 
	        neutral_advection_modes = [Ions.ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED], 
	        charged_diffusion_modes = [Ions.ION_CHARGED_DIFFUSION_MODE_PRESCRIBED], 
	        neutral_diffusion_modes = [Ions.ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED],
	        t_transp_expdecay_all_cs = t_exp_ion, t_transp_start_expdecay_all_cs = t_transp_onset, diffusion_initial_all_cs = DrrIon, advection_initial_all_cs = ArIon)
	else:
		ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardNe, Ninj=NinjNe, Zs=[10], isotopes=[0], molarFractions=[1], ionNames=['Ne_inj'], 
		abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxNe, shatterPoint=shatterPoint)


######################
# Run the simulation #
######################

# Include relevant data (except from the unknown quantities involved) in th output
ds.other.include('fluid', 'scalar')

# Name th output file
# Files without nShardD specified in the name use nShardD=300 as in the runs used for benchmarking with INDEX
# Files without NNe specified use NNe = 5e22
# Files without radius_wall specified use 2.15 m
filename_ending='ITER_'+scenario_name+'_INDEX_SPI_nShardD'+str(nShardD)+'NNe'+str(NinjD*molarFractionNe)+'abl_ioniz'+str(ds.eqsys.spi.abl_ioniz)+'radius_wall'+str(radius_wall)+'new'
if TQ_onset == 'Late':
    #filename_ending = filename_ending+'TQ_onset_from_output'
    filename_ending = filename_ending + 'late_TQ'
elif TQ_onset == 'Early':
    filename_ending = filename_ending + 'early_TQ'
if use_nre_transport or use_heat_transport:
    if transp_coeffs_filename is not None:
	    filename_ending = filename_ending + transp_coeffs_filename
    else:
	    filename_ending = filename_ending + '_t_duration'+str(t_duration)+'t_diffusion_over_t_duration'+str(t_duration_over_t_diffusion)
if use_ion_transport:
	filename_ending = filename_ending + 'iontransport_t_exp'+str(t_exp_ion)+'Ar'+str(ArIon)+'Drr'+str(DrrIon)
folder_name = 'output/'


# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)
	
if run_init:

	# Save settings to HDF5 file
	ds.save(folder_name+'init_settings'+filename_ending+'.h5')
	runiface(ds, folder_name+'output_init'+filename_ending+'.h5', quiet=False)


#######################
# RESTART set current #
#######################

# Extract neccessary flux surface averages from the output of the previous restart
# and use them to rescale the current density in the scenario data file
# to the current density at the outboard midplane, which should be given as input to DREAM.
# See presentation from ITER meeting 2022-03-07 for further details about this rescaling

do = DREAMOutput(folder_name+'output_init'+filename_ending+'.h5')

jprof=I.load_profile('J_Tor')
jprof = np.interp(do.grid.r, rref,jprof)
jprof = jprof * do.grid.GR0*do.grid.FSA_R02OverR2/do.grid.FSA_BOverBmin2/do.grid.Bmin/do.grid.R0
r_j = do.grid.r

I0=I.get_Ip() + I0_add_spike
ds.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_SELFCONSISTENT, inverse_wall_time = inverse_wall_time, R0=ds.radialgrid.R0)
ds.eqsys.j_ohm.setInitialProfile(j = jprof, radius = r_j, Ip0 = I0)


if run_init:
	# Save settings to HDF5 file
	ds.save(folder_name+'init_settings_'+filename_ending+'.h5')
	runiface(ds, folder_name+'output_init_set_current'+filename_ending+'.h5', quiet=False)

##########################
# RESTART injection init #
##########################
# Used to accelerate the distribution to carry the right current (only actually needed when using the kinetic solver)

ds2 = DREAMSettings(ds)
ds2.fromOutput(folder_name+'output_init_set_current'+filename_ending+'.h5')

ds2.timestep.setTmax(Tmax_init2)
ds2.timestep.setNt(Nt_init2)

if run_injection_init:
	ds2.save(folder_name+'init2_restart_settings_'+filename_ending+'.h5')
	runiface(ds2, folder_name+'output_init2_'+filename_ending+'.h5', quiet=False)



#####################
# RESTART injection #
#####################
# Here we actually make the first injection
ds3 = DREAMSettings(ds2)

# From now on, the temperature and electric field should be calculated self-consistently
ds3.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)
	
ds3.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds3.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_SELFCONSISTENT, inverse_wall_time = inverse_wall_time, R0=ds.radialgrid.R0)


ds3.fromOutput(folder_name+'output_init2_'+filename_ending+'.h5',ignore=['v_p'])

# Now make the shards from the first injection start moving
ds3.eqsys.spi.setShardVelocitiesUniform(nShard=None,abs_vp_mean=abs_vp_meanD,abs_vp_diff=abs_vp_diffD,alpha_max=alpha_maxD,nDim=2,add=False, shards=slice(0,nShardD))

# If the first injection contains any impurities, the pressure and current will be 
# significantly perturbed already during this injection stage, and then we
# turn on the transport when the shards reach the plasma boundary
if molarFractionNe>0:
	print('Turn on transport already during first injection')
	if use_heat_transport:
		T.setdBBFromDrr(ds3)
	
	if use_nre_transport:
		print(ds3.eqsys.n_re.transport.pstar)
		T.setSvenssonCoeff(ds3)
		if hotTailGrid_enabled:
		    Warning('Hot-tail grid enabled but no transport coefficients are set on the kinetic grid!')

ds3.timestep.setTmax(Tmax_injection)
ds3.timestep.setNt(Nt_injection)
ds3.timestep.setNumberOfSaveSteps(int(Tmax_injection/1e-4))

if run_injection:
	ds3.save(folder_name+'injection_restart_settings_'+filename_ending+'.h5')
	runiface(ds3, folder_name+'output_restart_injection_'+filename_ending+'.h5', quiet=False)
	
##############
# Restart CQ #
##############
ds4 = DREAMSettings(ds3)
ds4.fromOutput(folder_name+'output_restart_injection_'+filename_ending+'.h5',ignore=['x_p','v_p'])
ds4.timestep.setTmax(Tmax_CQ)
ds4.timestep.setNt(Nt_CQ)
ds4.timestep.setNumberOfSaveSteps(int(Tmax_CQ/1e-4))

# Turn on transport, if any
if use_heat_transport:
	T.t = T.t - ds3.timestep.tmax
	T.setdBBFromDrr(ds4)
	if not use_nre_transport:
	    T.shiftTimeSvensson(ds3, ds4)
if use_nre_transport:
	T.setSvenssonCoeff(ds4)
	if hotTailGrid_enabled:
	    Warning('Hot-tail grid enabled but no transport coefficients are set on the kinetic grid!')
	    
	T.shiftTimeSvensson(ds3, ds4)

if use_ion_transport:
	ds4.eqsys.n_i.ions[2].shiftTimeTranspCoeffs(ds3.timestep.tmax)
	ds4.eqsys.n_i.ions[3].shiftTimeTranspCoeffs(ds3.timestep.tmax)
	
# Set the shards of the second injection into motion and advance them until
# the fastest shards reach the plasma edge
if(nShardNe>0):
    do4=DREAMOutput(folder_name+'output_restart_injection_'+filename_ending+'.h5')
    ds4.eqsys.spi.vp=do4.eqsys.v_p.data[-1,:].flatten()
    ds4.eqsys.spi.xp=do4.eqsys.x_p.data[-1,:].flatten()
    ds4.eqsys.spi.setShardVelocitiesUniform(nShard=None, abs_vp_mean=abs_vp_meanNe, abs_vp_diff=abs_vp_diffNe, alpha_max=alpha_maxNe, nDim=2, add=False, shards=slice(-nShardNe,None))

    t_edge=(radius_wall-rref[-1])/np.max(-ds4.eqsys.spi.vp[-3*nShardNe::3])
    ds4.eqsys.spi.xp[-3*nShardNe:]=ds4.eqsys.spi.xp[-3*nShardNe:]+ds4.eqsys.spi.vp[-3*nShardNe:]*t_edge

if run_CQ:
	ds4.save(folder_name+'CQ_restart_settings_'+filename_ending+'.h5')
	runiface(ds4, folder_name+'output_restart_CQ_'+filename_ending+'.h5', quiet=False)
	
#################
# Restart CQ 2+ #
#################
ds5 = DREAMSettings(ds4)
for iCQ in range(nCQ_restart_start-2,nCQ_restart):
    if iCQ==0:
        ds5.fromOutput(folder_name+'output_restart_CQ_'+filename_ending+'.h5')
    else:
        ds5.fromOutput(folder_name+'output_restart_CQ'+str(iCQ+1)+'_'+filename_ending+'.h5')
    ds5.timestep.setTmax(Tmax_CQ_restart)
    ds5.timestep.setNt(Nt_CQ_restart)
    ds5.timestep.setNumberOfSaveSteps(int(Tmax_CQ_restart/1e-4))

    # Turn on transport, if any
    if use_heat_transport:
        T.t = T.t - ds4.timestep.tmax
        T.setdBBFromDrr(ds5)
        if not use_nre_transport:
            T.shiftTimeSvensson(ds4, ds5)
    if use_nre_transport:
        T.setSvenssonCoeff(ds5)
        if hotTailGrid_enabled:
            Warning('Hot-tail grid enabled but no transport coefficients are set on the kinetic grid!')

        T.shiftTimeSvensson(ds4, ds5)
		
    # Don't bother including the by now rather small ion transport coefficients, to decrease the runtime and need for time resolution 
    ds5.eqsys.n_i.ions[2].charged_diffusion_mode = Ions.ION_CHARGED_DIFFUSION_MODE_NONE
    ds5.eqsys.n_i.ions[2].neutral_diffusion_mode = Ions.ION_NEUTRAL_DIFFUSION_MODE_NONE
    ds5.eqsys.n_i.ions[2].charged_advection_mode = Ions.ION_CHARGED_ADVECTION_MODE_NONE
    ds5.eqsys.n_i.ions[2].neutral_advection_mode = Ions.ION_NEUTRAL_ADVECTION_MODE_NONE

    ds5.eqsys.n_i.ions[3].charged_diffusion_mode = Ions.ION_CHARGED_DIFFUSION_MODE_NONE
    ds5.eqsys.n_i.ions[3].neutral_diffusion_mode = Ions.ION_NEUTRAL_DIFFUSION_MODE_NONE
    ds5.eqsys.n_i.ions[3].charged_advection_mode = Ions.ION_CHARGED_ADVECTION_MODE_NONE
    ds5.eqsys.n_i.ions[3].neutral_advection_mode = Ions.ION_NEUTRAL_ADVECTION_MODE_NONE

    runiface(ds5, folder_name+'output_restart_CQ'+str(iCQ+2)+'_'+filename_ending+'.h5', quiet=False)
	

