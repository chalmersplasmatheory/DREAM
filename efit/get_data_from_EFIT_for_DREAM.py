''' 	The purpose of this script is to get data from EFIT
	specifically for SOFTv2 simulations
	and write to an hdf5 file
	See https://soft2.readthedocs.io/en/latest/magnetic_field.html#numerical-magnetic-field
	Author: Alex Tinguely 230817
	Usage: python3 get_data_from_EFIT_for_SOFT.py /path/to/gfile [/path/to/wallfile]
'''

### Import ###

import argparse
import numpy as np
from matplotlib import pyplot as plt
import scipy.interpolate as interp
from scipy.optimize import curve_fit
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from copy import *
import h5py
import sys
import geqdsk


### Setup ###

parser = argparse.ArgumentParser('Export data from EFIT to SOFT equilibrium file')

parser.add_argument('-H', '--no-hdf5', dest='hdf5',
                    help='Do NOT write EFIT data to HDF5 file.',
                    action='store_false', default=True)
parser.add_argument('-p', '--plot', help='Generate plots.',
                    action='store_true', default=False)
parser.add_argument('-v', '--verbose', help='Print to terminal.',
                    action='store_true', default=False)
parser.add_argument('-r', '--read-back',dest='read',
                    help='Try to read back the written file.',
                    action='store_true', default=False)
parser.add_argument('-t', '--itot',
                    help='Calculate total plasma current from Jphi.',
                    action='store_true', default=False)
parser.add_argument('-e', '--ienc',
                    help='Calculate enclosed plasma current from Br, Bz at LCFS.',
                    action='store_true', default=False)

parser.add_argument('gFilename', help='Name of EFIT Gfile to load equilibrium from.')
parser.add_argument('-w', '--wall', dest='wallfile',
                    help='Name of file containing wall.',
                    nargs='?', default='')

args = parser.parse_args()

### Define ###

#mf wall	Yes [1]		2-by-many vectorTokamak wall contour
def get_wall(wallfile):

	try:
		wall  	= np.loadtxt(wallfile,delimiter=',');
		if 	np.shape(wall)[0] != 2: wall = wall.T; 
		RFW  	= wall[0,:];       			# m, major radial points
		ZFW   	= wall[1,:];       		# m, vertical points
	except:
		raise ValueError('...Could not read user-supplied wall file');
		
	return wall,RFW,ZFW


def get_Br_Bz(psi,R,Z,dR=5e-3,dZ=5e-3,Rpsi=[],Zpsi=[]):

	''' 	Inputs
		- psi 	= poloidal flux [Wb]; matrix of size (R,Z)
		- R 	= array of major radii [m]; needs to be array or list
		- Z 	= array of vertical positions [m]; needs to be array or list
		- dR 	= major radial step size [m] to compute numerical derivative
		- dZ 	= vertical step size [m] to compute numerical derivative
		- Rpsi 	= major radial array for psi *if* not R
		- Zpsi 	= vertical array for psi *if* not Z
	'''

	# Check arguments
	if len(Rpsi)==0: Rpsi = R; 	# if not specified, set to R
	if len(Zpsi)==0: Zpsi = Z;	# if not specified, set to Z

	# Make positive (maybe not strictly necessary)
	dR 	= np.abs(dR);
	dZ 	= np.abs(dZ);

	# Initialize data
	BR 	= np.zeros(shape=(len(R),len(Z)));
	BZ 	= np.zeros(shape=(len(R),len(Z)));

	# Calculate (major) radial and vertical and toroidal B-fields
	''' 	In an axisymmetric toroidal geometry, we can relate Br and Bz to the flux psi in the following way:
		where here r is major radius and z is vertical position, usually R and Z

		Br(r,z) = -1/(2*pi*r)*Dpsi/Dz
		
		Bz(r,z) =  1/(2*pi*r)*Dpsi/Dr
	'''
	fpsi 	= interp.RectBivariateSpline(Rpsi,Zpsi,psi); # 2D spline function f_psi(R,Z)
	for iR in range(len(R)):
		for iZ in range(len(Z)):
			BR[iR,iZ] 	= -1/(2*np.pi*R[iR])*( fpsi(R[iR],	Z[iZ]+dZ/2)    	- fpsi(R[iR],	Z[iZ]-dZ/2)    	)/dZ; # T
			BZ[iR,iZ] 	=  1/(2*np.pi*R[iR])*( fpsi(R[iR]+dR/2,	Z[iZ]) 		- fpsi(R[iR]-dR/2,Z[iZ]   ) 	)/dR; # T
			
	return BR,BZ


### Get EFIT data ###

# Variable	Mandatory	Type		Description
#mf name	Yes		String	(Meta) 	Name of magnetic field data
#mf desc	Yes		String	(Meta) 	Description of data
name 		= args.gFilename;
desc 		= name;
gfile       	= geqdsk.geqdsk(name); 		# from Pablo RF


# Variable	Mandatory	Type		Description
#mf wall	Yes [1]		2-by-many vectorTokamak wall contour
if args.wallfile:	wall,RFW,ZFW 	= get_wall(args.wallfile);
else:
	RFW        	= gfile.Ginfo['RLIM'];       	# R of surrounding limiter contour in meter (FW=first wall)
	ZFW         	= gfile.Ginfo['ZLIM'];       	# Z of surrounding limiter contour in meter 
	wall 		= np.vstack((RFW,ZFW)); 	# wall points


#mf r		Yes				nr-vector	Radial grid
#mf z		Yes				nz-vector	Vertical grid
r       	= gfile.Ginfo['R'];		# m, major radial grid (array)
z       	= gfile.Ginfo['Z'];		# m, vertical position grid (array)		

RGRID,ZGRID 	= np.meshgrid(r,z); 		# make them 2D grids
dr		= np.mean(np.diff(r));  	# m, mean radial spacing
dz		= np.mean(np.diff(z));  	# m, mean vertical spacing


#mf Psi		No				nz-by-nr matrix	Poloidal magnetic flux
Psi        	= gfile.Ginfo['PSIRZ']*2*np.pi; # poloidal flux on 2D grid, psi(z,r) [Weber]; need factor of 2pi


#mf separatrix	Yes [1]				2-by-many vectorLast closed flux surface contour
RBBBS       	= gfile.Ginfo['RBBBS'];     	# m, major radial points at boundary
ZBBBS       	= gfile.Ginfo['ZBBBS'];     	# m, vertical points at boundary
separatrix	= np.vstack((RBBBS,ZBBBS)); 	# boundary of plasma

# Find points inside separatrix
BBBS 		= Polygon(np.vstack((RBBBS,ZBBBS)).T); 		# boundary of plasma
boolIn 		= np.zeros(shape=(len(r),len(z)),dtype=bool);
for ir in range(len(r)):
	for iz in range(len(z)):
		boolIn[ir,iz] = BBBS.contains(Point(r[ir],z[iz]));

'''
#mf Bphi	Yes				nz-by-nr matrix	Toroidal field component
B0		= gfile.Ginfo['BCENTR']; 	# Vacuum toroidal magnetic field in Tesla at RCENTR
R0	 	= gfile.Ginfo['RCENTR'];	# R in meter of vacuum toroidal magnetic field BCENTR
Bphi 		= B0*R0/RGRID; 			# T, toroidal magnetic field, Bphi(z,r)
BPHI		= Bphi.T; 			# Bphi(r,z)
'''

#mf Bphi	Yes				nz-by-nr matrix	Toroidal field component
FPOL		= gfile.Ginfo['FPOL']; 		# Poloidal current function in m-T, F = RBT on flux grid
Btor 		= FPOL/r; 			# T, toroidal magnetic field Btor(r)
Bphi 		= np.zeros(shape=np.shape(RGRID)); # to store data
for iZ in range(len(z)): Bphi[iZ,:] = Btor; 	# T, toroidal magnetic field, Bphi(z,r)
BPHI		= Bphi.T; 			# Bphi(r,z)


#mf Br		Yes				nz-by-nr matrix	Radial field component
#mf Bz		Yes				nz-by-nr matrix	Vertical field component
BR,BZ 		= get_Br_Bz(Psi.T,r,z); 	# T, Br(r,z); Bz(r,z)
Br 		= BR.T;				# T, Br(z,r)
Bz 		= BZ.T;				# T, Bz(z,r)

BR_IN 		= deepcopy(BR);			# T, Br(r,z)
BR_IN[~boolIn] 	= 0;				# T, Br(r,z), outside LCFS set to 0
BZ_IN 		= deepcopy(BZ);			# T, Bz(r,z)
BZ_IN[~boolIn] 	= 0;				# T, Bz(r,z), outside LCFS set to 0


#mf maxis	Yes				2-vector	Location of magnetic axis
RMAXIS       	= gfile.Ginfo['RMAXIS'];     	# m, major radius of magnetic axis
ZMAXIS       	= gfile.Ginfo['ZMAXIS'];     	# m, vertical position of magnetic axis
maxis		= np.vstack((RMAXIS,ZMAXIS));	# location of magnetic axis


#mf verBphi	No		nr-vector	Verification array for Bphi
#mf verBr	No		nr-vector	Verification array for Br
#mf verBz	No		nr-vector	Verification array for Bz
verBphi		= r;
verBr 		= r;
verBz 		= r;


if args.itot:# Get current density *inside LCFS*
	# From Ian Stewart email 230818, J_phi = R*p' + FF'/(R*mu_0)
	u0 		= 4*np.pi*1e-7; 
	PPRIME 		= gfile.Ginfo['PPRIMERZ']; 	# p'(z,r)
	FFPRIME 	= gfile.Ginfo['FFPRIMRZ'];	# FF'(z,r)
	Jphi 		= RGRID*PPRIME+FFPRIME/u0/RGRID;# A/m^2, jphi(z,r)
	JPHI		= Jphi.T; 			# A/m^2, jphi(r,z)
	JPHI_IN 	= deepcopy(JPHI);		# A/m^2, jphi(r,z)
	JPHI_IN[~boolIn]= 0; 				# A/m^2, jphi(r,z), outside LCFS set to 0
	ITOT 		= np.trapz(np.trapz(JPHI_IN,z,axis=1),r,axis=0); # integrate over z, then r
	if args.verbose: 	print('Itot = '+str(ITOT*1e-6)+' MA');


if args.ienc:# Integrate Bpol around LCFS, i.e. integral(B.dl) = u0*Ienc
	BR_LCFS 	= np.zeros(len(RBBBS)); 	# T, radial B-field at boundary
	BZ_LCFS 	= np.zeros(len(RBBBS));		# T, veritcal B-field at boundary
	for i in range(len(RBBBS)):
		BR_LCFS[i],BZ_LCFS[i] 	= get_Br_Bz(Psi.T,[RBBBS[i]],[ZBBBS[i]],Rpsi=r,Zpsi=z);
	integral1 	= np.sum(BR_LCFS[:-1]*np.diff(RBBBS)+BZ_LCFS[:-1]*np.diff(ZBBBS)); 	# calculate integral
	integral2 	= np.sum(BR_LCFS[1:]*np.diff(RBBBS)+BZ_LCFS[1:]*np.diff(ZBBBS));	# calculate integral again
	integral 	= 0.5*(integral1+integral2);						# take average
	Ienc 		= integral/u0; 			# A, enclosed current
	if args.verbose: 	print('Ienc = '+str(Ienc*1e-6)+' MA');


if args.hdf5:### Write to HDF5 ###
	
	h5fn 	= name+'_'*(args.wallfile!='')+args.wallfile+'.h5';
	with h5py.File(h5fn,'w') as hf:
		hf.create_dataset('Bphi',	data=Bphi); 	#	Yes		nz-by-nr matrix	Toroidal field component
		hf.create_dataset('Br',		data=Br); 	#	Yes		nz-by-nr matrix	Radial field component
		hf.create_dataset('Bz',		data=Bz); 	#	Yes		nz-by-nr matrix	Vertical field component
		hf.create_dataset('desc',	data=desc); 	#	Yes		String	(Meta) 	Description of data
		hf.create_dataset('maxis',	data=maxis); 	#	Yes		2-vector	Location of magnetic axis
		hf.create_dataset('name',	data=name); 	#	Yes		String	(Meta) 	Name of magnetic field data
		hf.create_dataset('Psi',	data=Psi); 	#	No		nz-by-nr matrix	Poloidal magnetic flux
		hf.create_dataset('r',		data=r); 	#	Yes		nr-vector	Radial grid
		hf.create_dataset('separatrix',	data=separatrix)#	Yes [1]		2-by-many vectorLast closed flux surface contour
		hf.create_dataset('verBphi',	data=verBphi); 	#	No		nr-vector	Verification array for Bphi
		hf.create_dataset('verBr',	data=verBr); 	#	No		nr-vector	Verification array for Br
		hf.create_dataset('verBz',	data=verBz); 	#	No		nr-vector	Verification array for Bz
		hf.create_dataset('wall',	data=wall); 	#	Yes [1]		2-by-many vectorTokamak wall contour
		hf.create_dataset('z',		data=z); 	#	Yes		nz-vector	Vertical grid
	
	if args.read:# Test getting file
		with h5py.File(h5fn,'r') as hfnew:
			keys 	= hfnew.keys();
			for key in keys:	print(key+', '+str(np.shape(hfnew.get(key))));
			
	
if args.plot:### Plots ###
	
	if 1: 	# Br, Bz
		fig 	= plt.figure();
		plt.contour(r,z,Psi,levels=21);
		plt.plot(RBBBS,ZBBBS,'r');
		plt.plot(RFW,ZFW,'k');
		plt.quiver(r,z,Br,Bz,pivot='mid',scale=30);
		#plt.colorbar();
		plt.xlabel('R (m)');
		plt.ylabel('Z (m)');
		plt.title('Radial and vertical B-field vectors');
		plt.axis('equal');
		fig.show();
		
	if 1: 	# Bphi
		fig 	= plt.figure();
		plt.contourf(r,z,Bphi);
		plt.plot(RBBBS,ZBBBS,'r');
		plt.plot(RFW,ZFW,'k');
		plt.colorbar();
		plt.xlabel('R (m)');
		plt.ylabel('Z (m)');
		plt.title('Toroidal B-field');
		plt.axis('equal');
		fig.show();
		
	if 1: 	# Psi
		fig 	= plt.figure();
		plt.contourf(r,z,Psi,levels=21);
		plt.plot(RBBBS,ZBBBS,'r');
		plt.plot(RFW,ZFW,'k');
		plt.colorbar();
		plt.xlabel('R (m)');
		plt.ylabel('Z (m)');
		plt.title('Poloidal flux');
		plt.axis('equal');
		fig.show();
		
	if 1: 	# JPHI
		fig 	= plt.figure();
		plt.contourf(r,z,Jphi,levels=21);
		plt.plot(RBBBS,ZBBBS,'r');
		plt.plot(RFW,ZFW,'k');
		plt.colorbar();
		plt.xlabel('R (m)');
		plt.ylabel('Z (m)');
		plt.title('Toroidal current density');
		plt.axis('equal');
		fig.show();

		
''' from https://soft2.readthedocs.io/en/latest/magnetic_field.html#numerical-magnetic-field
Variable	Mandatory	Type		Description
mf Bphi		Yes		nz-by-nr matrix	Toroidal field component
mf Br		Yes		nz-by-nr matrix	Radial field component
mf Bz		Yes		nz-by-nr matrix	Vertical field component
mf desc		Yes		String	(Meta) 	Description of data
mf maxis	Yes		2-vector	Location of magnetic axis
mf name		Yes		String	(Meta) 	Name of magnetic field data
mf Psi		No		nz-by-nr matrix	Poloidal magnetic flux
mf r		Yes		nr-vector	Radial grid
mf separatrix	Yes [1]		2-by-many vectorLast closed flux surface contour
mf verBphi	No		nr-vector	Verification array for Bphi
mf verBr	No		nr-vector	Verification array for Br
mf verBz	No		nr-vector	Verification array for Bz
mf wall		Yes [1]		2-by-many vectorTokamak wall contour
mf z		Yes		nz-vector	Vertical grid
[1]	(1, 2) At least one of the separatrix and wall variables must be present in the file.
'''

''' from G_EQDSK.pdf
RDIM: Horizontal dimension in meter of computational box
ZDIM: Vertical dimension in meter of computational box
RLEFT: Minimum R in meter of rectangular computational box
ZMID: Z of center of computational box in meter
RMAXIS: R of magnetic axis in meter
ZMAXIS: Z of magnetic axis in meter
SIMAG: poloidal flux at magnetic axis in Weber /rad
SIBRY: poloidal flux at the plasma boundary in Weber /rad
RCENTR: R in meter of vacuum toroidal magnetic field BCENTR
BCENTR: Vacuum toroidal magnetic field in Tesla at RCENTR
CURRENT: Plasma current in Ampere
FPOL: Poloidal current function in m-T, F = RBT on flux grid
PRES: Plasma pressure in nt / m2 on uniform flux grid
FFPRIM: FF'(y) in (mT)2
/ (Weber /rad) on uniform flux grid

PPRIME: P'(y) in (nt /m2

) / (Weber /rad) on uniform flux grid
PSIZR: Poloidal flux in Weber / rad on the rectangular grid points
QPSI: q values on uniform flux grid from axis to boundary
NBBBS: Number of boundary points
LIMITR: Number of limiter points
RBBBS: R of boundary points in meter
ZBBBS: Z of boundary points in meter
RLIM: R of surrounding limiter contour in meter
ZLIM: Z of surrounding limiter contour in meter
KVTOR Toroidal rotation switch
RVTOR Toroidal rotation characteristic major radius in m
NMASS Mass density switch
PRESSW Rotational pressure PW
in n/m2

PWPRIM PW

'(y) in (n/m2 )/(web/rad)

DMION Mass density on uniform poloidal flux grid
RHOVN Normalized toroidal flux on uniform poloidal flux grid
'''
