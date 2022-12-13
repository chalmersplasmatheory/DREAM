# This class provides routines for setting the electron transport coefficients (heat transport and runaway transport),
# either by reading the from a file or by estimating them with some of the available models. 
# Some functions in this class are modified versions of previous work done by Andréas Sundström

import numpy as np
import sys
import h5py as h5py
import os as os
import pathlib

path = str((pathlib.Path(__file__).parent / '../py/').resolve().absolute())
sys.path.append(path)

from DREAM.DREAMSettings import DREAMSettings
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMException import DREAMException
import DREAM.Settings.TransportSettings as Transport

## Elementary charge
c      = 299792458.      # m/s
e      = 1.602176634e-19 # C
m_e    = 9.1093837e-31   # kg
m_e_eV = 510.999e3       # eV


class TransportCoefficientReader:

    def __init__(self, filename = None, r_island = 0.0, coeff_scale = 1.0, t_before_onset = None, t_ramp = 1e-4, t_duration = None, Ip = None, interp3d=Transport.INTERP3D_LINEAR,interp1d=Transport.INTERP1D_LINEAR, interp1d_param=None):
        
        
        self.t = None # Time steps              [s]
        self.r = None # Radial coordinate       [m]
        self.xi = None # Pitch-angle coordinate  [1]
        self.p = None # Momentum coordinate     [m_e*c]
        self.Ar = None
        self.Drr = None
        
        self.Ip = Ip # Should be provided if one wants to interpolate the transport coefficients in plasma current instead of time
        
        self.interp3d = interp3d # Method for interpolating in r-xi-p
        self.interp1d = interp1d # Method for interpolating in t or Ip
        self.interp1d_param = interp1d_param # Which "time-variable" (Ip or t) to use
        
        if filename is not None:
            self.loadFile(filename = filename, r_island = r_island, coeff_scale = coeff_scale, t_before_onset = t_before_onset, t_ramp = t_ramp)
        
    def loadFile(self, filename, r_island = 0.0, coeff_scale = 1.0, t_before_onset = None, t_ramp = 1e-4):
        # Load transport coefficients and their coordinates, and amend according to input parameters
        
        with h5py.File(filename, 'r') as _f_coeff:
            self.t  = np.array(_f_coeff['time'])              # Time steps              [s]
            self.r  = np.array(_f_coeff['radius'])            # Radial coordinate       [m]
            self.xi = np.array(_f_coeff['pitch'])            # Pitch-angle coordinate  [1]
            self.p  = np.array(_f_coeff['momentum'])          # Momentum coordinate     [m_e*c]
            self.Ar  = np.array(_f_coeff['drift'])
            self.Drr = np.array(_f_coeff['diff'])  
            
        self.Ar  *= coeff_scale
        self.Drr *= coeff_scale
        
        # Set islands and time before onset, if any, and append zeros to turn of the transport after the desired duration    
        self.setRIsland(r_island)
        self.setTBeforeOnset(t_before_onset, t_ramp)
        self.appendZeros(t_ramp)
        
        
    def setCoeffsKinetic(self, dsObj, hotTailGridEnabled = True, runawayGridEnabled = False):
        # Set transport coefficients on the kinetic grid
        if hotTailGridEnabled:
            dsObj.eqsys.f_hot.transport.prescribedAdvection(self.Ar , t=self.t, r=self.r, p=self.p, xi=self.xi)
            dsObj.eqsys.f_hot.transport.prescribedDiffusion(self.Drr, t=self.t, r=self.r, p=self.p, xi=self.xi)
        if runawayGridEnabled:
            dsObj.eqsys.f_re.transport.prescribedAdvection(self.Ar , t=self.t, r=self.r, p=self.p, xi=self.xi)
            dsObj.eqsys.f_re.transport.prescribedDiffusion(self.Drr, t=self.t, r=self.r, p=self.p, xi=self.xi)

    def setSvenssonCoeff(self, dsObj):
        # Add the Svensson transport coefficients
        if self.interp1d_param is not None:
            if self.interp1d_param == Transport.SVENSSON_INTERP1D_PARAM_TIME and self.t is not None:
                dsObj.eqsys.n_re.transport.setSvenssonInterp1dParam(Transport.SVENSSON_INTERP1D_PARAM_TIME)
                dsObj.eqsys.n_re.transport.setSvenssonAdvection(self.Ar , t=self.t, r=self.r, p=self.p, xi=self.xi,
                                                                interp3d=self.interp3d, interp1d=self.interp1d)
                dsObj.eqsys.n_re.transport.setSvenssonDiffusion(self.Drr, t=self.t, r=self.r, p=self.p, xi=self.xi,
                                                                interp3d=self.interp3d, interp1d=self.interp1d)
                
            elif self.interp1d_param == Transport.SVENSSON_INTERP1D_PARAM_IP and self.Ip is not None:
                dsObj.eqsys.n_re.transport.setSvenssonInterp1dParam(Transport.SVENSSON_INTERP1D_PARAM_IP)
                dsObj.eqsys.n_re.transport.setSvenssonAdvection(self.Ar , Ip=self.Ip, r=self.r, p=self.p, xi=self.xi,
                                                                interp3d=self.interp3d, interp1d=self.interp1d)
                dsObj.eqsys.n_re.transport.setSvenssonDiffusion(self.Drr, Ip=self.Ip, r=self.r, p=self.p, xi=self.xi,
                                                                interp3d=self.interp3d, interp1d=self.interp1d)
            else:
                raise DREAMException("Error setting Svensson transport: interp1d_param is set but not its corresponding parameter.")

        elif self.Ip is not None and self.t is None:
            dsObj.eqsys.n_re.transport.setSvenssonInterp1dParam(Transport.SVENSSON_INTERP1D_PARAM_IP)
            dsObj.eqsys.n_re.transport.setSvenssonAdvection(self.Ar , Ip=self.Ip, r=self.r, p=self.p, xi=self.xi,
                                                            interp3d=self.interp3d, interp1d=self.interp1d)
            dsObj.eqsys.n_re.transport.setSvenssonDiffusion(self.Drr, Ip=self.Ip, r=self.r, p=self.p, xi=self.xi,
                                                            interp3d=self.interp3d, interp1d=self.interp1d)
        elif self.t is not None and self.Ip is None:
            dsObj.eqsys.n_re.transport.setSvenssonInterp1dParam(Transport.SVENSSON_INTERP1D_PARAM_TIME)
            dsObj.eqsys.n_re.transport.setSvenssonAdvection(self.Ar , t=self.t, r=self.r, p=self.p, xi=self.xi,
                                                            interp3d=self.interp3d, interp1d=self.interp1d)
            dsObj.eqsys.n_re.transport.setSvenssonDiffusion(self.Drr, t=self.t, r=self.r, p=self.p, xi=self.xi,
                                                            interp3d=self.interp3d, interp1d=self.interp1d)
        else:
            raise DREAMException("Error setting Svensson transport: the 1D-interpolation parameter has not been given.")
            
    def setdBBFromDrr(self, dsObj, ip = None, ixi = None, q = 1):
        # Calculates values of dB/B assuming Drr at momentum ndex ip 
        # and pitch angle index ixi corresponds to a Rechester-Rosenbluth transport coefficient,
        # and sets this dB/B for the heat transport of the settings object dsObj
        
        if ip is None:
            ip = np.argmin(self.p)
        if ixi is None:
            ixi = np.argmax(self.xi)
        dBB = np.sqrt(self.Drr[:,:,ixi,ip]/(np.pi*q*dsObj.radialgrid.R0*c*self.p[ip]/np.sqrt(1+self.p[ip]**2)*self.xi[ixi]))
        print(dBB)
        dsObj.eqsys.T_cold.transport.setMagneticPerturbation(dBB=dBB,r=self.r, t=self.t)
            
    
    def setRIsland(self, r_island = 0.0):
        # Helper function which helps putting in islands in the coeffs (or just set a zero transport coeff in the core if not specifiec)
        if self.r[0] != 0:
            if r_island > 0.0:
                self.r  = np.insert(self.r, 0, [0,self.r])
                _ins_array = np.zeros((2,self.t.size, self.xi.size, self.p.size))
            else:
                self.r  = np.insert(self.r, 0, [0])
                _ins_array = np.zeros((self.t.size, self.xi.size, self.p.size))
            # Adding in the leading zeros
            self.Ar  = np.insert(self.Ar,  0, _ins_array, axis=1)
            self.Drr = np.insert(self.Drr, 0, _ins_array, axis=1)
        else:
            if r_island > 0.0:
                self.Ar[:,self.r<r_island,:,:]  = 0.0
                self.Drr[:,self.r<r_island,:,:] = 0.0
                self.r  = np.insert(self.r, np.argmax(self.r>r_island), [r_island])
                _ins_array = np.zeros((self.t.size, self.xi.size, self.p.size))
                self.Ar  = np.insert(self.Ar,  0, _ins_array, axis=1)
                self.Drr = np.insert(self.Drr, 0, _ins_array, axis=1)
                
                
    def setTBeforeOnset(self, t_before_onset = None, t_ramp = 1e-4):
        # Helper function that adds zeros to the transport coefficients for times before the desired onset
        if t_before_onset is None:
            t_before_onset = self.t[0]
            
        self.t = self.t - self.t[0] + t_before_onset + t_ramp
            
        if t_before_onset>t_ramp:
            self.t  = np.insert(self.t, 0, [0,t_before_onset])
            _ins_array = np.zeros((2,self.r.size, self.xi.size, self.p.size))
        else:
            self.t  = np.insert(self.t, 0, [0])
            _ins_array = np.zeros((self.r.size, self.xi.size, self.p.size))
        # Adding in the leading zeros
        self.Ar  = np.insert(self.Ar,  0, _ins_array, axis=0)
        self.Drr = np.insert(self.Drr, 0, _ins_array, axis=0)
        if t_before_onset>0 and t_before_onset<t_ramp:
            print('WARNING: Time before onset of Svensson transport coefficients is specified but smaller than ramp-up time!')
            
        
    def appendZeros(self, t_ramp = 1e-4):
        # Helper function that appends two zeros in the time dimension to the transport coefficients, at times t_ramp and t_ramp+1 seconds
        # Can be used to turn off the transport at a given time, with a linear ramp down of duration t_ramp  , without stopping the simulation
           
        self.t  = np.append(self.t, [self.t[-1]+t_ramp, self.t[-1]+t_ramp+1])
        _ins_array = np.zeros((2,self.r.size, self.xi.size, self.p.size))

        # Adding in the appending zeros
        self.Ar  = np.append(self.Ar, _ins_array, axis=0)
        self.Drr = np.append(self.Drr, _ins_array, axis=0)
                
    
    def shiftTimeSvensson(self, ds_prev, ds_new):
        # Helper function that shifts the times of the coeffs
        ds_new.eqsys.n_re.transport.s_ar_t  = ds_prev.eqsys.n_re.transport.s_ar_t - ds_prev.timestep.tmax
        ds_new.eqsys.n_re.transport.s_drr_t = ds_prev.eqsys.n_re.transport.s_drr_t - ds_prev.timestep.tmax
        
    def tBeforeOnsetFromQCritAndPelletShardPosition(q, rref, shatterPoint, abs_vp_mean, qcrit = 2):
        # Estimates the time of the TQ onset, assuming the TQ starts when pellet shards 
        # starting from shatterPoint and traveling with velocity abs_vp_mean 
        # reach the flux surface where q=qcrit
        
        rqcrit = np.interp(qcrit, q, rref)
        return (shatterPoint[0] - rqcrit)/abs_vp_mean # NOTE: only valid for an injection close to Z=0!
    
    def tBeforeOnsetFromQCritAndTcoldFromOutput(q, qcrit, rref, Tcrit, filenames, folder):
        # Estimates the time of the TQ onset, assuming the TQ starts when the temperature
        # drops below Tcrit somewhere inside the flux surface where q=qcrit.
        # The temperature evolution used for determining this time i taken from the simulation output
        # named 'filename', which should be similar to the current simulation but without a transport-induced TQ
        
        if type(filenames) == str:
            filenames = [filenames]
            
        rcrit = np.interp(qcrit, q, rref)
        
        tBeforeCurrentRestart = 0
        for filename in filenames:
            doObj = DREAMOutput(folder+filename.strip())
            itTQHasStarted = np.any((doObj.eqsys.T_cold.data[:,:]<Tcrit)*(doObj.grid.r<rcrit),1)
            if np.any(itTQHasStarted):
                it = np.argwhere(itTQHasStarted)[0]
                return doObj.grid.t[it][0] + tBeforeCurrentRestart
            else:
                tBeforeCurentRestart = tBeforeCurrentRestart + doObj.grid.t[-1]
                
        raise DREAMException('Temperature does not drop below the critical temperature for the TQ onset inside the q=2 flux surface')
        
    def setDrrFromTQTimeScaleBesselApproximation(self, t_duration, a, Trepr, t_duration_over_t_diffusion = 1, t_before_onset = None, t_ramp = 1e-4, r_island = 0.0):
        # Estimate a radially flat diffusion coefficients based on a desired TQ duration t_duration. This is done assuming
        # that a fix ratio of the TQ duration and the diffusion time scale, t_duration_over_t_diffusion, gives
        # the desired temperature drop during the TQ duration. The diffusion time scale is calculated assuming a
        # bessel mode-like decay with the transport coefficient taken at the representative temperature Trepr.
        # The diffusion coefficient is then assumed to scale with momentum as p/(1+p^2)
        

        # Set up the grid on which the diffusion coefficients will be specified
        if t_duration is None:
            self.t = np.array([0])
        else:
            self.t = np.array([0,t_duration])
        self.r = np.array([0,a+1e-6])
        self.xi = np.array([-1,1])
        self.p = np.linspace(1e-6,100)
        
        # Assume no advection
        self.Ar = np.zeros((len(self.t), len(self.r), len(self.xi), len(self.p)))

        
        x1 = 2.4 # First zero of the bessel function (whichever it is...)
        vrepr = np.sqrt(2*e*Trepr/m_e) # thermal velocity corresponding to the representative temperature in m/s
        
        # Estimate the diffusion coefficient at the representative temperature
        Drr_vrepr = t_duration_over_t_diffusion * a**2/(t_duration * x1**2)
         
        # Set the p-dependent diffusion coefficient
        p_dep = self.p/(1+self.p**2)
        p_dep.reshape(1,1,1,-1)
        self.Drr = Drr_vrepr / vrepr * c * p_dep * np.ones((len(self.t), len(self.r), len(self.xi), len(self.p)))
        
        # Set islands and time before onset, if any, and append zeros to turn of the transport after the desired duration
        self.setRIsland(r_island)
        self.setTBeforeOnset(t_before_onset, t_ramp)
        self.appendZeros(t_ramp)
        
    def setPrecribedHyperResistivity(ds, Lambda0, Lambda1=0, t_before_onset = None, t_ramp = 1e-4, t_duration = None, t_sim_start = 0, r_island = 0.0, a_minor = None):
        if a_minor is None:
            a_minor = ds.radialgrid.a
        
        rLambda = np.array([0, a_minor])
        
        if t_before_onset is None:
            if t_duration is None:
                tLambda = np.array([0])
                Lambda_t_vec = np.array([Lambda0]).reshape(-1,1)
            else:
                tLambda = np.array([0, t_duration, t_ramp+t_duration, t_ramp+t_duration+1])
                Lambda_t_vec = np.array([Lambda0,Lambda0,Lambda1,Lambda1]).reshape(-1,1)
        else:
            if t_duration is None:
                tLambda = np.array([0,t_before_onset, t_before_onset+t_ramp, t_before_onset+t_ramp+1])
                Lambda_t_vec = np.array([0,0,Lambda0,Lambda0]).reshape(-1,1)   
            else:
                tLambda = np.array([0,t_before_onset, t_before_onset+t_ramp, t_before_onset+t_ramp+t_duration, t_before_onset+2*t_ramp+t_duration, t_before_onset+2*t_ramp+t_duration+1])
                Lambda_t_vec = np.array([0,0,Lambda0,Lambda0,Lambda1,Lambda1]).reshape(-1,1)        

        tLambda -= t_sim_start

        Lambda = Lambda_t_vec * np.ones((len(tLambda),len(rLambda)))
        ds.eqsys.psi_p.setHyperresistivity(Lambda, radius=rLambda, times=tLambda)
        
        

