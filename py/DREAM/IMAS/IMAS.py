# Functions to load data from IDS format and save DREAM output to IDS format

import imas
import numpy as np
from imas import imasdef
import logging
import os
import sys
#sys.path.append (os.environ['DREAMPATH']+'/py')
import h5py

import DREAM
from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
from scipy.interpolate import CubicSpline


def readInIDSSlice(shot, run, tokamak, user=os.getlogin(), time=-999, log=None, setUpDream=True, wall_radius=-999):
	
	# open shot file
	dataEntry = imas.DBEntry(imasdef.MDSPLUS_BACKEND, tokamak, int(shot), int(run), user_name=user)
	dataEntry.open()
	
	if log:
		logging.basicConfig(filename='%s.log' % log, format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG, filemode='w')
	
	# set time value if not specified by user
	if time==-999:
	
		try:
			eq_time = dataEntry.partial_get('equilibrium', 'time')
			cp_time = dataEntry.partial_get('core_profiles', 'time')
	
			time = max(min(eq_time), min(cp_time))
			if log:
				logging.info('The time used in loading IDS data is %s s\n', str(time))
		
		except:
			raise Exception('The equilibrium or the core profiles IDS seems to be empty. Please check the shot file!')
			
	# get the necessary IDS-s
	coreprof = dataEntry.get_slice('core_profiles', time, imasdef.LINEAR_INTERP)
	equilibrium = dataEntry.get_slice('equilibrium', time, imasdef.LINEAR_INTERP)
	wall = dataEntry.get('wall', 0)
	dataEntry.close()
	
	ion_names = []
	ion_densities = []
	ion_charges = []

	# load data from core profiles with error handling
	if len(coreprof.profiles_1d[0].grid.rho_tor) == 0:
		logging.warning('The minor radius from core profiles IDS seems to be empty!\n')
		logging.error('Core profiles IDS seems to be empty!\n')
	else:
		cp_radius = coreprof.profiles_1d[0].grid.rho_tor
		if log:
			logging.info('The minor radius from core profiles IDS loaded successfully!\n')

		if len(coreprof.profiles_1d[0].e_field_parallel) != len(cp_radius):
			logging.warning('The electric field from core profiles IDS seems to be empty!\n')
			e_field=-999
		else:
			e_field = coreprof.profiles_1d[0].e_field_parallel
			if log:
				logging.info('The electric field from core profiles IDS loaded successfully!\n')

		if len(coreprof.profiles_1d[0].electrons.temperature) != len(cp_radius):
			logging.warning('The electron temperature from core profiles IDS seems to be empty!\n')
			raise Exception('Error in loading core profiles data. Please check the input data!')
		else:
			T_cold = coreprof.profiles_1d[0].electrons.temperature
			if log:
				logging.info('The electron temperature from core profiles IDS loaded successfully!\n')

		if len(coreprof.profiles_1d[0].electrons.density) != len(cp_radius):
			logging.warning('The electron density from core profiles IDS seems to be empty!\n')
			raise Exception('Error in loading core profiles data. Please check the input data!')
		else:
			n_e = coreprof.profiles_1d[0].electrons.density
			if log:
				logging.info('The electron density from core profiles IDS loaded successfully!\n')

		for i in range(len(coreprof.profiles_1d[0].ion)):
			
			if len(coreprof.profiles_1d[0].ion[i].density) != len(cp_radius):
				logging.warning('An ion density from core profiles IDS seems to be empty!\n')
				raise Exception('Error in loading core profiles data. Please check the input data!')
			else:
				ion_names.append(coreprof.profiles_1d[0].ion[i].label)
				ion_densities.append(coreprof.profiles_1d[0].ion[i].density)
				ion_charges.append(int(coreprof.profiles_1d[0].ion[i].z_ion))
				
				if log:
					logging.info('An ion density from core profiles IDS loaded successfully!\n')
	
	if coreprof.vacuum_toroidal_field.r0 == -9e+40:
		logging.warning('The major radius of the magnetic axis from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		R = coreprof.vacuum_toroidal_field.r0
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')

	if len(coreprof.vacuum_toroidal_field.b0) == 0:
		logging.warning('The vacuum toroidal field from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		B0 = coreprof.vacuum_toroidal_field.b0
		
		if log:
			logging.info('The vacuum toroidal field from core profiles IDS loaded successfully!\n')
				
	#load equilibrium data with error handling
	if len(equilibrium.time_slice[0].profiles_1d.rho_tor) == 0:
		logging.warning('The minor radius from equilibrium IDS seems to be empty!\n')
		logging.error('Equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		eq_radius = equilibrium.time_slice[0].profiles_1d.rho_tor
		if log:
			logging.info('The minor radius from equilibrium IDS loaded successfully!\n')

		if len(equilibrium.time_slice[0].profiles_1d.j_tor) != len(eq_radius):
			logging.warning('The toroidal current density from equilibrium IDS seems to be empty!\n')
			raise Exception('Error in loading equilibrium data. Please check the input data!')
		else:
			j_tor = equilibrium.time_slice[0].profiles_1d.j_tor
			if log:
				logging.info('The toroidal current density from equilibrium IDS loaded successfully!\n')

	if equilibrium.time_slice[0].global_quantities.ip == 0:
		logging.warning('The plasma current from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		Ip = equilibrium.time_slice[0].global_quantities.ip
		if log:
			logging.info('The plasma current from equilibrium IDS loaded successfully!\n')
	
	# get the 2D psi data from equilibrium, and interpolate the core_profiles data onto this grid	
	# get the 2D grid
	if len(equilibrium.time_slice[0].profiles_2d[0].r) == 0:
		logging.warning('The 2D minor radius from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		r = equilibrium.time_slice[0].profiles_2d[0].r
		if log:
			logging.info('The 2D minor radius from equilibrium IDS loaded successfully!\n')

	if len(equilibrium.time_slice[0].profiles_2d[0].z) == 0:
		logging.warning('The 2D height coordinate from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		z = equilibrium.time_slice[0].profiles_2d[0].z
		if log:
			logging.info('The 2D height coordinate from equilibrium IDS loaded successfully!\n')
	
	# find the z = 0 column
	midplane_index = np.where(z==0)[1][0]
	
	# get the psi coordinates on the midplane
	if len(equilibrium.time_slice[0].profiles_2d[0].psi[:, midplane_index]) == 0:
		logging.warning('The 2D poloidal flux from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		eq_psi = equilibrium.time_slice[0].profiles_2d[0].psi[:, midplane_index]
		if log:
			logging.info('The 2D poloidal flux from equilibrium IDS loaded successfully!\n')	

	# limit the psi coordinates to the outer midplane
	# get the limiting psi values
	magnetic_axis_psi = max(eq_psi)			# psi value on the magnetic axis
	
	# get the indices fot the limits in psi_eq
	# get the magnetic axis index
	magnetic_axis_psi_index = np.where(eq_psi == max(eq_psi))[0][0]
	eq_psi = eq_psi[magnetic_axis_psi_index:]
	
	#normalize eq_psi
	eq_psi = eq_psi - eq_psi[0]	

	# load the coreprof psi
	if len(coreprof.profiles_1d[0].grid.psi) == 0:
		logging.warning('The poloidal flux from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		cp_psi = coreprof.profiles_1d[0].grid.psi
		if log:
			logging.info('The poloidal flux from core profiles IDS loaded successfully!\n')
	
	# interpolate the coreprof rho from the coreprof psi to the equilibrium psi
	interp_radius = CubicSpline(-cp_psi, cp_radius)
	radius = interp_radius(-eq_psi)
	
	# Check if the interpolation gave a correct answer
	error = (radius[-1] - cp_radius[-1]) / radius[-1]
	if abs(error) > 0.1:
		raise Exception("The interpolation gave significant difference in the equilibrium and coreprofile radii. Cannot continue, please check the input parameters.")
	
	#load wall data with error handling
	if wall.description_2d[0].vessel.unit[0].annular.resistivity == -9e+40:
		logging.warning('The the wall resistivity from wall IDS seems to be empty!\n')
		raise Exception('Error in loading wall data. Please check the input data!')
	else:
		resistivity = wall.description_2d[0].vessel.unit[0].annular.resistivity
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')

	if setUpDream:
		
		ds = DREAMSettings()
		
		ds.radialgrid.setCustomGridPoints(radius)

		if wall_radius == -999:
			raise Exception('If you want to set up a dream settings object, please specify the wall radius!')
		else:
			ds.radialgrid.setWallRadius(wall_radius)
			
		ds.radialgrid.setMajorRadius(R)
		ds.radialgrid.setB0(abs(B0))
		
		for i in range(len(coreprof.profiles_1d[0].ion)):
			ds.eqsys.n_i.addIon(ion_names[i], Z=ion_charges[i], iontype=Ions.IONS_DYNAMIC, Z0=ion_charges[i], n=ion_densities[i], r=cp_radius)
			
		ds.eqsys.T_cold.setInitialProfile(T_cold, radius=cp_radius)
		ds.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
		
		if e_field!=-999:
			ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
			ds.eqsys.E_field.setInitialProfile(efield=e_field, radius=cp_radius)
		
		else:
			inverse_wall_time = resistivity / (4*np.pi*1e-7*R*np.log(R/wall_radius))
			if log:
				logging.info('The inverse wall time is %s 1/s with the given wall radius.\n', inverse_wall_time)
		
			ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
			ds.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_SELFCONSISTENT, inverse_wall_time=inverse_wall_time, R0=R)
		
		ds.eqsys.j_ohm.setInitialProfile(j_tor, radius=eq_radius, Ip0=Ip)
		
		return ds
	
	else:
		
		IDS_variables = {'radius': radius, 'Core_profiles_radius': cp_radius, 'Equilibrium_radius': eq_radius, 'wall_radius': wall_radius, 'Major_radius': R, 'Magnetic_field': B0, 'Ion_names': ion_names, 'Ion_charges': ion_charges, 'Ion_densities': ion_densities, 'Temperature': T_cold, 'Electric_field': e_field,  'Wall_resistivity': resistivity, 'Toroidal_j': j_tor, 'Plasma_current': Ip, 'Time': time}

		return IDS_variables
	
	
def readInIDS(shot, run, tokamak, user=os.getlogin(), log=False, setUpDream=True, wall_radius=-999):
	
	#full IDS read in
	
	# open shot file
	dataEntry = imas.DBEntry(imasdef.MDSPLUS_BACKEND, tokamak, int(shot), int(run), user_name=user)
	dataEntry.open()
	
	if log:
		logging.basicConfig(filename='IDS_summary.log', format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG, filemode='w')
		
	# get the necessary IDS-s
	#coreprof = dataEntry.get_slice('core_profiles')
	#equilibrium = dataEntry.get_slice('equilibrium')
	#wall = dataEntry.get('wall',0)
	
	ion_names = []
	ion_densities = []
	ion_charges = []
	
	# load data from core profiles with error handling
	# coreprof radius
	cp_radius = dataEntry.partial_get('core_profiles', 'profiles_1d(:)/grid/rho_tor')
	
	if cp_radius.ndim == 1:
		logging.warning('The minor radius from core profiles IDS seems to be empty!\n')
		logging.error('Core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		if log:
			logging.info('The minor radius from core profiles IDS loaded successfully!\n')

	# coreprof electric field
	cp_efield = dataEntry.partial_get('core_profiles', 'profiles_1d(:)/e_field_parallel')
	
	if cp_efield.ndim == 1:
		logging.warning('The electric field from core profiles IDS seems to be empty!\n')
		e_field=-999
	else:
		if log:
			logging.info('The electric field from core profiles IDS loaded successfully!\n')

	# coreprof electron temperature
	cp_Te = np.transpose(dataEntry.partial_get('core_profiles', 'profiles_1d(:)/electrons/temperature'))
	
	if cp_Te.ndim == 1:
		logging.warning('The electron temperature from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		if log:
			logging.info('The electron temperature from core profiles IDS loaded successfully!\n')

				
	# coreprof electron density
	cp_ne = dataEntry.partial_get('core_profiles', 'profiles_1d(:)/electrons/density')
	
	if cp_ne.ndim == 1:
		logging.warning('The electron density from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		if log:
			logging.info('The electron density from core profiles IDS loaded successfully!\n')			

	# coreprof ion densities
	cp_ions = dataEntry.partial_get('core_profiles', 'profiles_1d(:)/ion')
	
	if cp_ions.ndim == 0:
		logging.warning('There are no ion species in the core profiles IDS!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		if log:
			logging.info('The ions species from core profiles IDS loaded successfully!\n')	
	
	for i in range(np.shape(cp_ions)[0]):
		
		# get the ion names and charges
		ion_names.append(cp_ions[i,0].label)
		ion_charges.append(int(cp_ions[i,0].z_ion))
				
		for j in range(np.shape(cp_ions)[1]):
						
			if len(cp_ions[i,j].density) == 0:
				logging.warning('An ion density from core profiles IDS seems to be empty!\n')
				raise Exception('Error in loading core profiles data. Please check the input data!')
			else:
				ion_densities.append(cp_ions[i,j].density)
				
				if log:
					logging.info('An ion density from core profiles IDS loaded successfully!\n')
	
	# coreprof major radius
	R = dataEntry.partial_get('core_profiles', 'vacuum_toroidal_field/r0')
	
	if R == -9e+40:
		logging.warning('The major radius of the magnetic axis from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')

	# coreprof toroidal field
	B0 = dataEntry.partial_get('core_profiles', 'vacuum_toroidal_field/b0')
	
	if len(B0) == 0:
		logging.warning('The vacuum toroidal field from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		if log:
			logging.info('The vacuum toroidal field from core profiles IDS loaded successfully!\n')

	# coreprof time
	cp_time = dataEntry.partial_get('core_profiles', 'time')
	
	if len(cp_time) == 0:
		logging.warning('The time from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		if log:
			logging.info('The time from core profiles IDS loaded successfully!\n')

	#load equilibrium data with error handling
	# equilibrium radius
	eq_radius = dataEntry.partial_get('equilibrium', 'time_slice(:)/profiles_1d/rho_tor')
	
	if eq_radius.ndim == 1:
		logging.warning('The minor radius from equilibrium IDS seems to be empty!\n')
		logging.error('Equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		if log:
			logging.info('The minor radius from equilibrium IDS loaded successfully!\n')
	
	# equilibrium toroidal current density
	eq_jtor = dataEntry.partial_get('equilibrium', 'time_slice(:)/profiles_1d/j_tor')

	if eq_jtor.ndim == 1:
		logging.warning('The toroidal current density from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		if log:
			logging.info('The toroidal current density from equilibrium IDS loaded successfully!\n')

	# equilibrium total current
	eq_ip = dataEntry.partial_get('equilibrium', 'time_slice(:)/global_quantities/ip')
	
	if len(eq_ip) == 0:
		logging.warning('The plasma current from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		if log:
			logging.info('The plasma current from equilibrium IDS loaded successfully!\n')
			
	# equilibrium time
	eq_time = dataEntry.partial_get('equilibrium', 'time')
	
	if len(eq_time) == 0:
		logging.warning('The time from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibirum data. Please check the input data!')
	else:
		if log:
			logging.info('The time from equilibrium IDS loaded successfully!\n')
	
	# get the first timepoint where both IDSs are filled
	time = max(min(eq_time), min(cp_time))
	if log:
		logging.info('The time used in loading IDS data is %s s\n', str(time))
	
	
	# get the 2D psi data from equilibrium, and interpolate the core_profiles data onto this grid
	# since DREAM cannot handle time varying radial grids, the interpolation will be done on the first time where both IDS-s are available, and keep that radial grid for the full simulation
	# get the 2D grid
	# get the necessary IDSs
	coreprof = dataEntry.get_slice('core_profiles', time, imasdef.LINEAR_INTERP)
	equilibrium = dataEntry.get_slice('equilibrium', time, imasdef.LINEAR_INTERP)
	wall = dataEntry.get('wall', 0)
	
	if len(equilibrium.time_slice[0].profiles_2d[0].r) == 0:
		logging.warning('The 2D minor radius from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		r = equilibrium.time_slice[0].profiles_2d[0].r
		if log:
			logging.info('The 2D minor radius from equilibrium IDS loaded successfully!\n')

	if len(equilibrium.time_slice[0].profiles_2d[0].z) == 0:
		logging.warning('The 2D height coordinate from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		z = equilibrium.time_slice[0].profiles_2d[0].z
		if log:
			logging.info('The 2D height coordinate from equilibrium IDS loaded successfully!\n')
	
	# find the z = 0 column
	midplane_index = np.where(z==0)[1][0]
	
	# get the psi coordinates on the midplane
	if len(equilibrium.time_slice[0].profiles_2d[0].psi[:, midplane_index]) == 0:
		logging.warning('The 2D poloidal flux from equilibrium IDS seems to be empty!\n')
		raise Exception('Error in loading equilibrium data. Please check the input data!')
	else:
		eq_psi = equilibrium.time_slice[0].profiles_2d[0].psi[:, midplane_index]
		if log:
			logging.info('The 2D poloidal flux from equilibrium IDS loaded successfully!\n')	

	# limit the psi coordinates to the outer midplane
	# get the limiting psi values
	magnetic_axis_psi = max(eq_psi)			# psi value on the magnetic axis
	
	# get the indices fot the limits in psi_eq
	# get the magnetic axis index
	magnetic_axis_psi_index = np.where(eq_psi == max(eq_psi))[0][0]
	eq_psi = eq_psi[magnetic_axis_psi_index:]
	
	# normalize eq_psi
	eq_psi = eq_psi - eq_psi[0]	

	# load the coreprof psi
	if len(coreprof.profiles_1d[0].grid.psi) == 0:
		logging.warning('The poloidal flux from core profiles IDS seems to be empty!\n')
		raise Exception('Error in loading core profiles data. Please check the input data!')
	else:
		cp_psi = coreprof.profiles_1d[0].grid.psi
		if log:
			logging.info('The poloidal flux from core profiles IDS loaded successfully!\n')
	
	# load the coreprof radius
	single_cp_radius = coreprof.profiles_1d[0].grid.rho_tor
	
	# interpolate the coreprof rho from the coreprof psi to the equilibrium psi
	interp_radius = CubicSpline(-cp_psi, single_cp_radius)
	radius = interp_radius(-eq_psi)
	
	# Check if the interpolation gave a correct answer
	error = (radius[-1] - single_cp_radius[-1]) / radius[-1]
	if abs(error) > 0.1:
		raise Exception("The interpolation gave significant difference in the equilibrium and coreprofile radii. Cannot continue, please check the input parameters.")
	
	# load wall data with error handling
	if wall.description_2d[0].vessel.unit[0].annular.resistivity == -9e+40:
		logging.warning('The the wall resistivity from wall IDS seems to be empty!\n')
		raise Exception('Error in loading wall data. Please check the input data!')
	else:
		resistivity = wall.description_2d[0].vessel.unit[0].annular.resistivity
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')

	dataEntry.close()

	if setUpDream:
		
		ds = DREAMSettings()
		
		ds.radialgrid.setCustomGridPoints(radius)

		if wall_radius == -999:
			raise Exception('If you want to set up a dream settings object, please specify the wall radius!')
		else:
			ds.radialgrid.setWallRadius(wall_radius)
			
		ds.radialgrid.setMajorRadius(R)
		ds.radialgrid.setB0(abs(B0[0]))
		
		cp_time_length = len(cp_time)
		
		for i in range(len(ion_names)):
				
			ds.eqsys.n_i.addIon(ion_names[i], Z=ion_charges[i], iontype=Ions.IONS_PRESCRIBED, Z0=ion_charges[i], n=ion_densities[i*cp_time_length:i*cp_time_length + cp_time_length], r=cp_radius[:,0], t=cp_time)
			
		ds.eqsys.T_cold.setPrescribedData(cp_Te, radius=cp_radius[:,0], times=cp_time)
		ds.eqsys.T_cold.setType(Tcold.TYPE_PRESCRIBED)
		
		if e_field!=-999:
			ds.eqsys.E_field.setType(ElectricField.TYPE_PRESCRIBED)
			ds.eqsys.E_field.setPrescribedData(efield=cp_efield, radius=cp_radius[:,0], times=cp_time)
		
		else:
			inverse_wall_time = resistivity / (4*np.pi*1e-7*R*np.log(R/wall_radius))
			if log:
				logging.info('The inverse wall time is %s 1/s with the given wall radius.\n', inverse_wall_time)
		
			ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
			ds.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_SELFCONSISTENT, inverse_wall_time=inverse_wall_time, R0=R)
		
		ds.eqsys.j_ohm.setInitialProfile(eq_jtor[:,0], radius=eq_radius[:,0], Ip0=eq_ip[0])
		
		return ds
	
	else:
		IDS_variables = {'radius': radius, 'Core_profiles_radius': cp_radius, 'Equilibrium_radius': eq_radius, 'wall_radius': wall_radius, 'Major_radius': R, 'Magnetic_field': B0, 'Ion_names': ion_names, 'Ion_charges': ion_charges, 'Ion_densities': ion_densities, 'Temperature': T_cold, 'Electric_field': e_field,  'Wall_resistivity': resistivity, 'Toroidal_j': j_tor, 'Plasma_current': Ip, 'Time': time}

		return IDS_variables
