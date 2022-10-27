# Load data from IDS format and save DREAM output to IDS format

import imas
import numpy as np
from imas import imasdef
import logging
import os


def readInIDS(shot, run, tokamak, user=os.getlogin(), time=-999, log=False, setUpDream=True):
	
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

	if wall.description_2d[0].vessel.unit[0].annular.resistivity == -9e+40:
		logging.warning('The the wall resistivity from wall IDS seems to be empty!\n')
	else:
		resistivity = wall.description_2d[0].vessel.unit[0].annular.resistivity
		if log:
			logging.info('The major radius of the magnetic axis from core profiles IDS loaded successfully!\n')

	return 'DONE'

