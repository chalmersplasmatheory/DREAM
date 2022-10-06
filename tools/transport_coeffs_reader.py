# This class provides routines for setting the electron transport coefficients (heat transport and runaway transport),
# either by reading the from a file or by estimating them with some of the available models. 
# Some functions in this class are modified versions of previous work done by Andreas Sundström

import numpy as np
import sys
import h5py as h5py
import os as os

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMException import DREAMException
import DREAM.Settings.TransportSettings as Transport

## Elementary charge
c      = 299792458.      # m/s
e      = 1.602176634e-19 # C
m_e    = 9.1093837e-31   # kg
m_e_eV = 510.999e3       # eV


class transport_coeffs_reader:

    ####### Helper functions for setting the Svensson coeffs.#############################
    def __init__(self, filename = None, r_island = 0.0, coeff_scale = 1.0, t_before_onset = None, t_ramp = 1e-4, t_duration = None, Ip = None, interp3d=Transport.INTERP3D_LINEAR,interp1d=Transport.INTERP1D_LINEAR, interp1d_param=None):
        ############### Load transport coefficients and their coordinates, and amend according to input parameters ##########
        
        if filename is None:
            # If no filename is specified, assume that there is no advection
            # but a homogenious diffusion coefficient with an assumed p-dependence of p/(1+p^2).
            # By setting an arbitrary diffusion coefficient with this scaling here,
            # the shift in time of the transport onset is properly taken care of later in this constructor,
            # and then the magnitude of the diffusion coefficient can simply be rescaled later 
            # (see for example the setDrrFromTQTimeScaleBesselApproximation-function)
            
            if t_duration is None:
                t_coeff = np.array([0])
            else:
                t_coeff = np.array([0,t_duration])
            r_coeff = np.array([0,10])
            xi_coeff = np.array([-1,1])
            p_coeff = np.linspace(1e-6,100)
            
            p_dep = p_coeff/(1+p_coeff**2)
            p_dep.reshape(1,1,1,-1)
            
            Ar_coeff = np.zeros((len(t_coeff), len(r_coeff), len(xi_coeff), len(p_coeff)))
            Drr_coeff = np.ones((len(t_coeff), len(r_coeff), len(xi_coeff), len(p_coeff)))*p_dep
        else:
            try:
                _f_coeff  = h5py.File(filename, 'r')
            except:
                _f_coeff  = h5py.File(os.path.dirname(__file__)+'/'+filename, 'r')

            t_coeff  = np.array(_f_coeff['time'])              # Time steps              [s]
            r_coeff  = np.array(_f_coeff['radius'])            # Radial coordinate       [m]
            xi_coeff = np.array(_f_coeff['pitch'])            # Pitch-angle coordinate  [1]
            p_coeff  = np.array(_f_coeff['momentum'])          # Momentum coordinate     [m_e*c]
            Ar_coeff  = np.array(_f_coeff['drift'])
            Drr_coeff = np.array(_f_coeff['diff'])  

            _f_coeff.close()

        r_coeff, Ar_coeff, Drr_coeff = transport_coeffs_reader.rIslandSvenssonCoeff(t_coeff, r_coeff, xi_coeff, p_coeff, Ar_coeff, Drr_coeff, r_island)
        t_coeff, Ar_coeff, Drr_coeff = transport_coeffs_reader.tBeforeOnsetSvenssonCoeff(t_coeff, r_coeff, xi_coeff, p_coeff, Ar_coeff, Drr_coeff, t_before_onset, t_ramp)
        t_coeff, Ar_coeff, Drr_coeff = transport_coeffs_reader.appendZerosSvenssonCoeff(t_coeff, r_coeff, xi_coeff, p_coeff, Ar_coeff, Drr_coeff, t_ramp)
        
        Ar_coeff  *= coeff_scale
        Drr_coeff *= coeff_scale
        
        self.t = t_coeff
        self.r = r_coeff
        self.xi = xi_coeff
        self.p = p_coeff
        self.Ar = Ar_coeff
        self.Drr = Drr_coeff
        
        self.Ip = Ip
        
        self.interp3d = interp3d
        self.interp1d = interp1d
        self.interp1d_param = interp1d_param
        
    def setDrr(self, dsObj, hotTailGridEnabled = True, runawayGridEnabled = False):
        if hotTailGridEnabled:
            dsObj.eqsys.f_hot.transport.prescribedAdvection(self.Ar , t=self.t, r=self.r, p=self.p, xi=self.xi)
            dsObj.eqsys.f_hot.transport.prescribedDiffusion(self.Drr, t=self.t, r=self.r, p=self.p, xi=self.xi)
        if runawayGridEnabled:
            dsObj.eqsys.f_re.transport.prescribedAdvection(self.Ar , t=self.t, r=self.r, p=self.p, xi=self.xi)
            dsObj.eqsys.f_re.transport.prescribedDiffusion(self.Drr, t=self.t, r=self.r, p=self.p, xi=self.xi)

    def setSvenssonCoeff(self, dsObj):
        # Add the Svensson transport coefficients
        if self.interp1d_param is not None:
            if self.interp1d_param == Transport.SVENSSON_INTERP1D_PARAM_TIME and self.t_coeff is not None:
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
            
    
    def rIslandSvenssonCoeff(t_coeff, r_coeff, xi_coeff, p_coeff, Ar_coeff, Drr_coeff, r_island = 0):
        # Helper function which helps putting in islands in the coeffs (or just set a zero transport coeff in the core if not specifiec)
        if r_coeff[0] != 0:
            if r_island > 0.0:
                r_coeff  = np.insert(r_coeff, 0, [0,r_island])
                _ins_array = np.zeros((2,t_coeff.size, xi_coeff.size, p_coeff.size))
            else:
                r_coeff  = np.insert(r_coeff, 0, [0])
                _ins_array = np.zeros((t_coeff.size, xi_coeff.size, p_coeff.size))
            # Adding in the leading zeros
            Ar_coeff  = np.insert(Ar_coeff,  0, _ins_array, axis=1)
            Drr_coeff = np.insert(Drr_coeff, 0, _ins_array, axis=1)
        else:
            if r_island > 0.0:
                Ar_coeff[:,r_coeff<r_island,:,:]  = 0.0
                Drr_coeff[:,r_coeff<r_island,:,:] = 0.0
                r_coeff  = np.insert(r_coeff, np.argmax(r_coeff>r_island), [r_island])
                _ins_array = np.zeros((t_coeff.size, xi_coeff.size, p_coeff.size))
                Ar_coeff  = np.insert(Ar_coeff,  0, _ins_array, axis=1)
                Drr_coeff = np.insert(Drr_coeff, 0, _ins_array, axis=1)
                
        return r_coeff, Ar_coeff, Drr_coeff
                
    def tBeforeOnsetSvenssonCoeff(t_coeff, r_coeff, xi_coeff, p_coeff, Ar_coeff, Drr_coeff, t_before_onset = None, t_ramp = 1e-4):
        # Adds zeros to the transport coefficients for times before the desired onset
        if t_before_onset is None:
            t_before_onset = t_coeff[0]
            
        t_coeff = t_coeff - t_coeff[0] + t_before_onset + t_ramp
            
        if t_before_onset>t_ramp:
            t_coeff  = np.insert(t_coeff, 0, [0,t_before_onset])
            _ins_array = np.zeros((2,r_coeff.size, xi_coeff.size, p_coeff.size))
        else:
            t_coeff  = np.insert(t_coeff, 0, [0])
            _ins_array = np.zeros((r_coeff.size, xi_coeff.size, p_coeff.size))
        # Adding in the leading zeros
        Ar_coeff  = np.insert(Ar_coeff,  0, _ins_array, axis=0)
        Drr_coeff = np.insert(Drr_coeff, 0, _ins_array, axis=0)
        if t_before_onset>0 and t_before_onset<t_ramp:
            print('WARNING: Time before onset of Svensson transport coefficients is specified but smaller than ramp-up time!')
            
        return t_coeff, Ar_coeff, Drr_coeff
        
    def appendZerosSvenssonCoeff(t_coeff, r_coeff, xi_coeff, p_coeff, Ar_coeff, Drr_coeff, t_ramp = 1e-4):
            
        t_coeff  = np.append(t_coeff, [t_coeff[-1]+t_ramp, t_coeff[-1]+t_ramp+1])
        _ins_array = np.zeros((2,r_coeff.size, xi_coeff.size, p_coeff.size))

        # Adding in the appending zeros
        Ar_coeff  = np.append(Ar_coeff, _ins_array, axis=0)
        Drr_coeff = np.append(Drr_coeff, _ins_array, axis=0)
            
        return t_coeff, Ar_coeff, Drr_coeff
                
    
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
                return doObj.grid.t[it] + tBeforeCurrentRestart
            else:
                tBeforeCurentRestart = tBeforeCurrentRestart + doObj.grid.t[-1]
                
        raise DREAMException('Temperature does not drop below the critical temperature for the TQ onset inside the q=2 flux surface')
        
    def setDrrFromTQTimeScaleBesselApproximation(self, t_duration, a, Trepr, t_duration_over_t_diffusion = 1):
        # Estimate the diffusion coefficients based on a desired TQ duration t_duration. This is done assuming
        # that a fix value of the TQ duration and the diffusion time scale, t_duration_over_t_diffusion, gives
        # the desired temperature drop during the TQ duration. The diffusion time scale is calculated assuming a
        # bessel mode-like decay with the transport coefficient taken at the representative temperature Trepr.
        # The diffusion coefficient is then assumed to scale with momentum as p/(1+p^2)
        
        x1 = 2.4 # First zero of the bessel function (whichever it is...)
        vrepr = np.sqrt(2*e*Trepr/m_e) # thermal velocity corresponding to the representative temperature in m/s
        
        # Estimate the diffusion coefficient at the representative temperature
        Drr_vrepr = t_duration_over_t_diffusion * a**2/(t_duration * x1**2)
        
        # Find time indices where Drr>0 (set in the constructor)
        it = np.argwhere(self.Drr[:,0,0,0]>0)[0] 
        
        # Rescale the diffusion coefficients>0 to give the desired value at the representative temperature
        self.Drr *= Drr_vrepr / vrepr * c * self.p[0]/(1+self.p[0]**2) / self.Drr[it,0,0,0]
        
        
        

