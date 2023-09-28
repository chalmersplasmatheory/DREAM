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


def readInIDSSlice(shot, run, tokamak, user=os.getlogin(), time=-999, log=False, setUpDream=True, wall_radius=-999):
	
	# open shot file
	dataEntry = imas.DBEntry(imasdef.MDSPLUS_BACKEND, tokamak, int(shot), int(run), user_name=user)
	dataEntry.open()
	
	if log:
		logging.basicConfig(filename='IDS_summary.log', format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG, filemode='w')
	
	# set time value if not specified by user
	if time==-999:
	
		try:
			timeEq = dataEntry.partial_get('equilibrium', 'time')
			timeCp = dataEntry.partial_get('core_profiles', 'time')
	
			time = max(min(timeEq), min(timeCp))
			if log:
				logging.info('The time used in loading IDS data is %s s\n', str(time))
		
		except:
			raise Exception('The equilibrium or the core profiles IDS seems to be empty. Please check the shot file!')
			
	# get the necessary IDS-s
	coreprof = dataEntry.get_slice('core_profiles', time, imasdef.LINEAR_INTERP)
	equilibrium = dataEntry.get_slice('equilibrium', time, imasdef.LINEAR_INTERP)
	wall = dataEntry.get('wall',0)
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
		else:
			T_cold = coreprof.profiles_1d[0].electrons.temperature
			if log:
				logging.info('The electron temperature from core profiles IDS loaded successfully!\n')

		if len(coreprof.profiles_1d[0].electrons.density) != len(cp_radius):
			logging.warning('The electron density from core profiles IDS seems to be empty!\n')
		else:
			n_e = coreprof.profiles_1d[0].electrons.density
			if log:
				logging.info('The electron density from core profiles IDS loaded successfully!\n')

		for i in range(len(coreprof.profiles_1d[0].ion)):
			
			if len(coreprof.profiles_1d[0].ion[i].density) != len(cp_radius):
				logging.warning('An ion density from core profiles IDS seems to be empty!\n')
			else:
				ion_names.append(coreprof.profiles_1d[0].ion[i].label)
				ion_densities.append(coreprof.profiles_1d[0].ion[i].density)
				ion_charges.append(int(coreprof.profiles_1d[0].ion[i].z_ion))
				
				if log:
					logging.info('An ion density from core profiles IDS loaded successfully!\n')
	
	if coreprof.vacuum_toroidal_field.r0 == -9e+40:
		logging.warning('The major radius of the magnetic axis from core profiles IDS seems to be empty!\n')
	else:
		R = coreprof.vacuum_toroidal_field.r0
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')

		if len(coreprof.vacuum_toroidal_field.b0) == 0:
			logging.warning('The vacuum toroidal field from core profiles IDS seems to be empty!\n')
		else:
			B0 = coreprof.vacuum_toroidal_field.b0
			if log:
				logging.info('The vacuum toroidal field from core profiles IDS loaded successfully!\n')
				
	#load equilibrium data with error handling
	if len(equilibrium.time_slice[0].profiles_1d.rho_tor) == 0:
		logging.warning('The minor radius from equilibrium IDS seems to be empty!\n')
		logging.error('Equilibrium IDS seems to be empty!\n')
	else:
		eq_radius = equilibrium.time_slice[0].profiles_1d.rho_tor
		if log:
			logging.info('The minor radius from equilibrium IDS loaded successfully!\n')

		if len(equilibrium.time_slice[0].profiles_1d.j_tor) != len(eq_radius):
			logging.warning('The toroidal current density from equilibrium IDS seems to be empty!\n')
		else:
			j_tor = equilibrium.time_slice[0].profiles_1d.j_tor
			if log:
				logging.info('The toroidal current density from equilibrium IDS loaded successfully!\n')

	if equilibrium.time_slice[0].global_quantities.ip == 0:
		logging.warning('The plasma current from equilibrium IDS seems to be empty!\n')
	else:
		Ip = equilibrium.time_slice[0].global_quantities.ip
		if log:
			logging.info('The plasma current from equilibrium IDS loaded successfully!\n')

	#load wall data with error handling
	if wall.description_2d[0].vessel.unit[0].annular.resistivity == -9e+40:
		logging.warning('The the wall resistivity from wall IDS seems to be empty!\n')
	else:
		resistivity = wall.description_2d[0].vessel.unit[0].annular.resistivity
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')
	
	# get the 2D psi data from equilibrium, and interpolate the core_profiles data onto this grid
	
	# get the 2D grid
	r = equilibrium.time_slice[0].profiles_2d[0].r
	z = equilibrium.time_slice[0].profiles_2d[0].z
	
	# find the z = 0 column
	midplane_index = np.where(z==0)[1][0]
	
	# get the psi coordinates on the midplane
	eq_psi = equilibrium.time_slice[0].profiles_2d[0].psi[:, midplane_index]
	
	# limit the psi coordinates to the outer midplane
	
	# get the limiting psi values
	magnetic_axis_psi = max(eq_psi)			# psi value on the magnetic axis
	
	# get the indices fot the limits in psi_eq
	# get the magnetic axis index
	magnetic_axis_psi_index = np.where(eq_psi == max(eq_psi))[0][0]
	eq_psi = eq_psi[magnetic_axis_psi_index:]
	
	#normalize eq_psi
	eq_psi = eq_psi - eq_psi[0	

	# interpolate the coreprof rho from the coreprof psi to the equilibrium psi
	cp_psi = coreprof.profiles_1d[0].grid.psi

	interp_radius = CubicSpline(-cp_psi, cp_radius)
	radius = interp_radius(-eq_psi)
	
	# Check if the interpolation gave a correct answer
	error = (radius[-1] - cp_radius[-1]) / radius[-1]
	if error > 0.1:
		raise Exception("The interpolation gave significant difference in the equilibrium and coreprofile radii. Cannot continue, please check the input parameters.")

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
				logging.info('The iverse wall time is %s 1/s with the given wall radius.\n', inverse_wall_time)
		
			ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
			ds.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_SELFCONSISTENT, inverse_wall_time=inverse_wall_time, R0=R)
		
		ds.eqsys.j_ohm.setInitialProfile(j_tor, radius=eq_radius, Ip0=Ip)
		
		return ds
	
	else:
		print('return a list? hdf5?')
		return 'TO DO'
	
	
def readInIDS(shot, run, tokamak, user=os.getlogin(), log=False, setUpDream=True, wall_radius=-999):
	
	# TODO full IDS read in
	
	# open shot file
	dataEntry = imas.DBEntry(imasdef.MDSPLUS_BACKEND, tokamak, int(shot), int(run), user_name=user)
	dataEntry.open()
	
	if log:
		logging.basicConfig(filename='IDS_summary.log', format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG, filemode='w')
		
	# get the necessary IDS-s
	coreprof = dataEntry.get_slice('core_profiles')
	equilibrium = dataEntry.get_slice('equilibrium')
	wall = dataEntry.get('wall',0)
	dataEntry.close()
	
	ion_names = []
	ion_densities = []
	ion_charges = []
	
	# load data from core profiles with error handling
	cp_radius = dataEntry.partial_get('core_profiles', 'profiles_1d[(0)/grid/rho_tor')
	
	if len(coreprof.profiles_1d[0].grid.rho_tor) == 0:
		logging.warning('The minor radius from core profiles IDS seems to be empty!\n')
		logging.error('Core profiles IDS seems to be empty!\n')
		
	elif len(coreprof.profiles_1d[0].grid.rho_tor) == 0:
		logging.warning('The minor radius from core profiles IDS seems to be empty!\n')
		logging.error('Core profiles IDS seems to be empty!\n')
	else:

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
		else:
			T_cold = coreprof.profiles_1d[0].electrons.temperature
			if log:
				logging.info('The electron temperature from core profiles IDS loaded successfully!\n')

		if len(coreprof.profiles_1d[0].electrons.density) != len(cp_radius):
			logging.warning('The electron density from core profiles IDS seems to be empty!\n')
		else:
			n_e = coreprof.profiles_1d[0].electrons.density
			if log:
				logging.info('The electron density from core profiles IDS loaded successfully!\n')

		for i in range(len(coreprof.profiles_1d[0].ion)):
			
			if len(coreprof.profiles_1d[0].ion[i].density) != len(cp_radius):
				logging.warning('An ion density from core profiles IDS seems to be empty!\n')
			else:
				ion_names.append(coreprof.profiles_1d[0].ion[i].label)
				ion_densities.append(coreprof.profiles_1d[0].ion[i].density)
				ion_charges.append(int(coreprof.profiles_1d[0].ion[i].z_ion))
				
				if log:
					logging.info('An ion density from core profiles IDS loaded successfully!\n')
	
	if coreprof.vacuum_toroidal_field.r0 == -9e+40:
		logging.warning('The major radius of the magnetic axis from core profiles IDS seems to be empty!\n')
	else:
		R = coreprof.vacuum_toroidal_field.r0
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')

		if len(coreprof.vacuum_toroidal_field.b0) == 0:
			logging.warning('The vacuum toroidal field from core profiles IDS seems to be empty!\n')
		else:
			B0 = coreprof.vacuum_toroidal_field.b0
			if log:
				logging.info('The vacuum toroidal field from core profiles IDS loaded successfully!\n')
	
