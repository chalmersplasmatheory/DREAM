
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from numpy.matlib import repmat
from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.IonSpecies import IonSpecies, IONS_PRESCRIBED, IONIZATION_MODE_FLUID, IONIZATION_MODE_KINETIC, IONIZATION_MODE_KINETIC_APPROX_JAC, ION_OPACITY_MODE_TRANSPARENT, ION_CHARGED_DIFFUSION_MODE_NONE, ION_CHARGED_DIFFUSION_MODE_PRESCRIBED, ION_NEUTRAL_DIFFUSION_MODE_NONE, ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED, ION_CHARGED_ADVECTION_MODE_NONE, ION_CHARGED_ADVECTION_MODE_PRESCRIBED, ION_NEUTRAL_ADVECTION_MODE_NONE, ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED
from . UnknownQuantity import UnknownQuantity
from .. import AdvectionInterpolation

# Model to use for ion heat
IONS_T_I_NEGLECT = 1
IONS_T_I_INCLUDE = 2

# Interpolation methods for advection term in transport equation
AD_INTERP_CENTRED  = AdvectionInterpolation.AD_INTERP_CENTRED
AD_INTERP_UPWIND   = AdvectionInterpolation.AD_INTERP_UPWIND
AD_INTERP_UPWIND_2ND_ORDER = AdvectionInterpolation.AD_INTERP_UPWIND_2ND_ORDER
AD_INTERP_DOWNWIND = AdvectionInterpolation.AD_INTERP_DOWNWIND
AD_INTERP_QUICK    = AdvectionInterpolation.AD_INTERP_QUICK 
AD_INTERP_SMART    = AdvectionInterpolation.AD_INTERP_SMART 
AD_INTERP_MUSCL    = AdvectionInterpolation.AD_INTERP_MUSCL 
AD_INTERP_OSPRE    = AdvectionInterpolation.AD_INTERP_OSPRE 
AD_INTERP_TCDF     = AdvectionInterpolation.AD_INTERP_TCDF  

AD_INTERP_JACOBIAN_LINEAR = AdvectionInterpolation.AD_INTERP_JACOBIAN_LINEAR
AD_INTERP_JACOBIAN_FULL   = AdvectionInterpolation.AD_INTERP_JACOBIAN_FULL  
AD_INTERP_JACOBIAN_UPWIND = AdvectionInterpolation.AD_INTERP_JACOBIAN_UPWIND

class Ions(UnknownQuantity):
    

    def __init__(self, settings, ionization=IONIZATION_MODE_FLUID):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.ions = list()
        self.r    = None
        self.t    = None
        
        self.rChargedPrescribedDiffusion = None
        self.rNeutralPrescribedDiffusion = None
        self.tChargedPrescribedDiffusion = None
        self.tNeutralPrescribedDiffusion = None
        
        self.rChargedPrescribedAdvection = None
        self.rNeutralPrescribedAdvection = None
        self.tChargedPrescribedAdvection = None
        self.tNeutralPrescribedAdvection = None
        self.tSourceTerm                 = None

        self.ionization = ionization
        self.typeTi = IONS_T_I_NEGLECT
        
        self.advectionInterpolationCharged = AdvectionInterpolation.AdvectionInterpolation(kinetic=False)
        self.advectionInterpolationNeutral = AdvectionInterpolation.AdvectionInterpolation(kinetic=False)


    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, Z0=None, isotope=0, SPIMolarFraction=-1, opacity_mode=ION_OPACITY_MODE_TRANSPARENT, 
        charged_diffusion_mode=ION_CHARGED_DIFFUSION_MODE_NONE, charged_prescribed_diffusion=None, rChargedPrescribedDiffusion=None, tChargedPrescribedDiffusion=None,
        neutral_diffusion_mode=ION_NEUTRAL_DIFFUSION_MODE_NONE, neutral_prescribed_diffusion=None, rNeutralPrescribedDiffusion=None, tNeutralPrescribedDiffusion=None,
        charged_advection_mode=ION_CHARGED_ADVECTION_MODE_NONE, charged_prescribed_advection=None, rChargedPrescribedAdvection=None, tChargedPrescribedAdvection=None,
        neutral_advection_mode=ION_NEUTRAL_ADVECTION_MODE_NONE, neutral_prescribed_advection=None, rNeutralPrescribedAdvection=None, tNeutralPrescribedAdvection=None,
        t_transp_expdecay_all_cs = None, t_transp_start_expdecay_all_cs = 0, diffusion_initial_all_cs = None, diffusion_final_all_cs = 0, advection_initial_all_cs = None, advection_final_all_cs = 0, r_expdecay_all_cs = None, t_expdecay_all_cs = None, 
        init_equil=False, T=None, n=None, r=None, t=None, tritium=False, hydrogen=False):

        """
        Adds a new ion species to the plasma.

        :param str name:        Name by which the ion species will be referred to.
        :param int Z:           Ion charge number.
        :param int isotope:            Ion mass number.
        :param int iontype:     Method to use for evolving ions in time.
        :param int Z0:          Charge state to populate (used for populating exactly one charge state for the ion).
        :param n:               Ion density (can be either a scalar, 1D array or 2D array, depending on the other input parameters)
        :param float SPIMolarFraction: Molar fraction of the SPI injection (if any). A negative value means that this species is not part of the SPI injection 
        :param numpy.ndarray r: Radial grid on which the input density is defined.
        :param T:               Ion initial temperature (can be scalar for uniform temperature, otherwise 1D array matching `r` in size)
        :param numpy.ndarray r: Radial grid on which the input density and temperature is defined.
        :param numpy.ndarray t: Time grid on which the input density is defined.
        :param bool tritium:    If ``True``, the ion species is treated as Tritium.
        :param bool hydrogen:   If ``True``, the ion species is treated as Hydrogen (single proton).
        """
        if (self.r is not None) and (r is not None) and (np.any(self.r[:] != r[:])):
            if self.r.size == 1:
                self.changeRadialGrid(r)
            else:
                raise EquationException("The radial grid must be the same for all ion species.")
        if (self.t is not None) and (t is not None) and (np.any(self.t != t)):
            raise EquationException("The time grid must be the same for all ion species.")
            
        if (self.rChargedPrescribedDiffusion is not None) and (rChargedPrescribedDiffusion is not None) and (np.any(self.rChargedPrescribedDiffusion != rChargedPrescribedDiffusion)):
            raise EquationException("The radial grid for the prescribed charged diffusion must be the same for all ion species.")
        if (self.tChargedPrescribedDiffusion is not None) and (tChargedPrescribedDiffusion is not None) and (np.any(self.tChargedPrescribedDiffusion != tChargedPrescribedDiffusion)):
            raise EquationException("The time grid for the prescribed charged diffusion must be the same for all ion species.")
            
        if (self.rNeutralPrescribedDiffusion is not None) and (rNeutralPrescribedDiffusion is not None) and (np.any(self.rNeutralPrescribedDiffusion != rNeutralPrescribedDiffusion)):
            raise EquationException("The radial grid for the prescribed neutral diffusion must be the same for all ion species.")
        if (self.tNeutralPrescribedDiffusion is not None) and (tNeutralPrescribedDiffusion is not None) and (np.any(self.tNeutralPrescribedDiffusion != tNeutralPrescribedDiffusion)):
            raise EquationException("The time grid for the prescribed neutral diffusion must be the same for all ion species.")

        if T is not None:
            self.typeTi = IONS_T_I_INCLUDE

        ion = IonSpecies(settings=self.settings, name=name, Z=Z, ttype=iontype, Z0=Z0, isotope=isotope, SPIMolarFraction=SPIMolarFraction, opacity_mode=opacity_mode, 
            charged_diffusion_mode=charged_diffusion_mode, charged_prescribed_diffusion=charged_prescribed_diffusion, rChargedPrescribedDiffusion=rChargedPrescribedDiffusion, tChargedPrescribedDiffusion=tChargedPrescribedDiffusion,
            neutral_diffusion_mode=neutral_diffusion_mode, neutral_prescribed_diffusion=neutral_prescribed_diffusion, rNeutralPrescribedDiffusion=rNeutralPrescribedDiffusion, tNeutralPrescribedDiffusion=tNeutralPrescribedDiffusion,           
            charged_advection_mode=charged_advection_mode, charged_prescribed_advection=charged_prescribed_advection, rChargedPrescribedAdvection=rChargedPrescribedAdvection, tChargedPrescribedAdvection=tChargedPrescribedAdvection,
            neutral_advection_mode=neutral_advection_mode, neutral_prescribed_advection=neutral_prescribed_advection, rNeutralPrescribedAdvection=rNeutralPrescribedAdvection, tNeutralPrescribedAdvection=tNeutralPrescribedAdvection,
            t_transp_expdecay_all_cs = t_transp_expdecay_all_cs, t_transp_start_expdecay_all_cs = t_transp_start_expdecay_all_cs,
            diffusion_initial_all_cs = diffusion_initial_all_cs, diffusion_final_all_cs = diffusion_final_all_cs, 
            advection_initial_all_cs = advection_initial_all_cs, advection_final_all_cs = advection_final_all_cs, 
            r_expdecay_all_cs = r_expdecay_all_cs, t_expdecay_all_cs = t_expdecay_all_cs,            
            init_equil=init_equil, T=T, n=n, r=r, t=t, interpr=self.r, interpt=None, tritium=tritium, hydrogen=hydrogen)

        self.ions.append(ion)

        self.r = ion.getR()
        if ion.getTime() is not None:
            self.t = ion.getTime()
            
        if charged_diffusion_mode==ION_CHARGED_DIFFUSION_MODE_PRESCRIBED:
            self.rChargedPrescribedDiffusion = ion.getRChargedPrescribedDiffusion()
            self.tChargedPrescribedDiffusion = ion.getTChargedPrescribedDiffusion()
            
        if neutral_diffusion_mode==ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED:
            self.rNeutralPrescribedDiffusion = ion.getRNeutralPrescribedDiffusion()
            self.tNeutralPrescribedDiffusion = ion.getTNeutralPrescribedDiffusion()
            
        if charged_advection_mode==ION_CHARGED_ADVECTION_MODE_PRESCRIBED:
            self.rChargedPrescribedAdvection = ion.getRChargedPrescribedAdvection()
            self.tChargedPrescribedAdvection = ion.getTChargedPrescribedAdvection()
            
        if neutral_advection_mode==ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED:
            self.rNeutralPrescribedAdvection = ion.getRNeutralPrescribedAdvection()
            self.tNeutralPrescribedAdvection = ion.getTNeutralPrescribedAdvection()


    def addIonSource(self, species, dNdt=None, t=None, Z0=0):
        """
        Add a source term for the specified ion species.

        :param species: Name of the ion species to add this source term to.
        :param dNdt:    Number of particles to add per unit time.
        :param t:       Time grid associated with ``dNdt`` (if any).
        :param Z0:      For scalar or 1D ``dNdt``, the charge state for which to add the source term.
        """
        if t is None:
            t = self.tSourceTerm
        elif self.tSourceTerm is not None and not np.all(t == self.tSourceTerm):
            raise EquationException(f"The time grid used for ion sources must be the same for all ion species.")
            
        found = False
        for ion in self.ions:
            if ion.name == species:
                ion.initialize_source(n=dNdt, t=t, Z0=Z0)
                found = True

                self.tSourceTerm = ion.getSourceTime()

        if not found:
            raise EquationException(f"No ion species with name '{species}' has been added to the simulation. Unable to add source term.")


    def changeRadialGrid(self, r):
        """
        Change the radial grid used for the ion species.
        """
        for ion in self.ions:
            if ion.r.size == 1:
                ion.n = ion.n * np.ones(ion.n.shape[:-1] + (r.size,))
                ion.T = ion.T * np.ones(ion.T.shape[:-1] + (r.size,))
                ion.r = r
            else:
                fn = scipy.interpolate.interp1d(ion.r, ion.n, axis=-1, bounds_error=False, fill_value='extrapolate')
                fT = scipy.interpolate.interp1d(ion.r, ion.n, axis=-1, bounds_error=False, fill_value='extrapolate')
                ion.n = fn(r)
                ion.T = fT(r)
                ion.r = r

        self.r = r


    def getCharges(self):
        """
        Returns a list of the charges of the various ion species
        contained by this object.
        """
        return [ion.getZ() for ion in self.ions]


    def getIsotopes(self):
        """
        Returns a list of the isotopes of the various ion species
        contained by this object.
        """
        return [ion.getIsotope() for ion in self.ions]


    def getSPIMolarFraction(self):
        """
        Returns a list of the SPI molar fractions of the various ion species
        contained by this object.
        """
        return [ion.getSPIMolarFraction() for ion in self.ions]


    def getIndex(self, species):
        """
        Return the index of the ion species with the specified name.
        """
        for iIon in range(len(self.ions)):
            if self.ions[iIon].getName() == species:
                return iIon

        raise EquationException(f"No species with name '{species}'.")


    def getIon(self, i=None):
        """
        Returns the ion species with the specified index or name.

        :param i: Index or name of ion species to retrieve.
        """
        if type(i) == int: return self.ions[i]
        elif type(i) == str:
            for j in range(0, len(self.ions)):
                if self.ions[j].getName() == i:
                    return self.ions[j]

            raise EquationException("No ion with name '{}' has been defined.".format(i))
        else:
            raise EquationException("Invalid call to 'getIon()'.")


    def setIonization(self, ionization=IONIZATION_MODE_FLUID):
        """
        Sets which model to use for ionization.

        :param int ionization: Flag indicating which model to use for ionization.
        """
        self.ionization=ionization


    def getHydrogenSpecies(self):
        """
        Returns a list of names of the ion species which are treated
        as Hydrogen.
        """
        hydr = []
        for ion in self.ions:
            if ion.hydrogen:
                hydr.append(ion.getName())

        return hydr


    def getTritiumSpecies(self):
        """
        Returns a list of names of the ion species which are treated
        as Tritium.
        """
        trit = []
        for ion in self.ions:
            if ion.tritium:
                trit.append(ion.getName())

        return trit


    def getTypes(self):
        """
        Returns a list of ion types for the various ion species
        contained by this object.
        """
        return [ion.getType() for ion in self.ions]
        
    def getOpacityModes(self):
        """
        Returns a list of ion opacity modes for the various ion species
        contained by this object.
        """
        return [ion.getOpacityMode() for ion in self.ions]
        
    def getChargedDiffusionModes(self):
        """
        Returns a list of ion charged diffusion modes for the various ion species
        contained by this object.
        """
        return [ion.getChargedDiffusionMode() for ion in self.ions]
        
    def getNeutralDiffusionModes(self):
        """
        Returns a list of ion neutral diffusion modes for the various ion species
        contained by this object.
        """
        return [ion.getNeutralDiffusionMode() for ion in self.ions]
        
    def getChargedAdvectionModes(self):
        """
        Returns a list of ion charged advection modes for the various ion species
        contained by this object.
        """
        return [ion.getChargedAdvectionMode() for ion in self.ions]
        
    def getNeutralAdvectionModes(self):
        """
        Returns a list of ion neutral advection modes for the various ion species
        contained by this object.
        """
        return [ion.getNeutralAdvectionMode() for ion in self.ions]
        
    def setAdvectionInterpolationMethodCharged(self, ad_int=AD_INTERP_CENTRED,
        ad_jac=AD_INTERP_JACOBIAN_FULL, fluxlimiterdamping=1.0):
        """
        Sets the interpolation method that is used in the charged advection terms of
        the transport equation.
        
        :param int ad_int:               Interpolation method to use for the radial coordinate.
        :param int ad_jac:               Jacobian interpolation mode to use for the radial coordinate.
        :param float fluxlimiterdamping: Damping parameter used to under-relax the interpolation coefficients during non-linear iterations (should be between 0 and 1).
        """
        self.advectionInterpolationCharged.setMethod(ad_int=ad_int, ad_jac=ad_jac, fluxlimiterdamping=fluxlimiterdamping)
        
    def setAdvectionInterpolationMethodNeutral(self, ad_int=AD_INTERP_CENTRED,
        ad_jac=AD_INTERP_JACOBIAN_FULL, fluxlimiterdamping=1.0):
        """
        Sets the interpolation method that is used in the neutral advection terms of
        the transport equation.
        
        :param int ad_int:               Interpolation method to use for the radial coordinate.
        :param int ad_jac:               Jacobian interpolation mode to use for the radial coordinate.
        :param float fluxlimiterdamping: Damping parameter used to under-relax the interpolation coefficients during non-linear iterations (should be between 0 and 1).
        """
        self.advectionInterpolationNeutral.setMethod(ad_int=ad_int, ad_jac=ad_jac, fluxlimiterdamping=fluxlimiterdamping)


    def setIonType(self, index, ttype):
        """
        Modifies the type of equation used for the specified ion species.

        :param index: Index or name of ion species to set type for.
        :param int ttype: Type of equation to use for evolving the ion species.
        """
        ion = self.getIon(index)

        # Note that the DREAM kernel only uses positive type indices.
        # The negative type indices are interface extensions which can
        # only be used with the 'initialize()' methods.
        if ttype <= 0:
            raise DREAMException("Trying to set invalid ion type for ion species '{}': {}.".format(ion.name, ttype))

        ion.ttype = ttype

        
    def shiftTimeTranspCoeffs(self, tShift):
        """
        Shift the time grids for the ion transport coefficients by an amount tShift. This is needed between restarts.
        
        :param float tShift: Amount of time the time grids for the ion transport coefficient should be shifted
        """
        if self.tChargedPrescribedDiffusion is not None:
            self.tChargedPrescribedDiffusion = self.tChargedPrescribedDiffusion - tShift
        if self.tNeutralPrescribedDiffusion is not None:    
            self.tNeutralPrescribedDiffusion = self.tNeutralPrescribedDiffusion - tShift
        if self.tChargedPrescribedAdvection is not None:
            self.tChargedPrescribedAdvection = self.tChargedPrescribedAdvection - tShift
        if self.tNeutralPrescribedAdvection is not None:
            self.tNeutralPrescribedAdvection = self.tNeutralPrescribedAdvection - tShift
    

    def fromdict(self, data):
        """
        Load settings from the specified dictionary.
        
        :param dict data: Dictionary containing all settings to load.
        """
        names        = data['names'].split(';')[:-1]
        Z            = data['Z']

        if 'isotopes' in data and len(data['isotopes']) == len(Z):
            isotopes = data['isotopes']
        else:
            isotopes = [0]*len(Z)

        if 'types' in data and len(data['types']) == len(Z):
            types = data['types']
        else:
            types = [0]*len(Z)

        if 'opacity_modes' in data and len(data['opacity_modes']) == len(Z):
            opacity_modes = data['opacity_modes']
        else:
            opacity_modes = self.getOpacityModes()

        charged_diffusion_modes = [ION_CHARGED_DIFFUSION_MODE_NONE]*len(Z)
        neutral_diffusion_modes = [ION_NEUTRAL_DIFFUSION_MODE_NONE]*len(Z)
        charged_advection_modes = [ION_CHARGED_ADVECTION_MODE_NONE]*len(Z)
        neutral_advection_modes = [ION_NEUTRAL_ADVECTION_MODE_NONE]*len(Z)

        if 'charged_diffusion_modes' in data:
            charged_diffusion_modes = data['charged_diffusion_modes']
        if 'neutral_diffusion_modes' in data:
            neutral_diffusion_modes = data['neutral_diffusion_modes']
        if 'charged_advection_modes' in data:
            charged_advection_modes = data['charged_advection_modes']
        if 'neutral_advection_modes' in data:
            neutral_advection_modes = data['neutral_advection_modes']

        SPIMolarFraction = data['SPIMolarFraction']
        nZSPI = len(Z)-np.sum(SPIMolarFraction<0)
        if nZSPI>0:
            nShard = int(np.sum(SPIMolarFraction>=0)/nZSPI)
        else:
            nShard = 0

        if 'tritiumnames' in data:
            tritiumnames = data['tritiumnames'].split(';')[:-1]
        else:
            tritiumnames = []

        if 'hydrogennames' in data:
            hydrogennames = data['hydrogennames'].split(';')[:-1]
        else:
            hydrogennames = []

        initial    = None
        prescribed = None
        initialTi  = None
        charged_prescribed_diffusion = None
        neutral_prescribed_diffusion = None
        charged_prescribed_advection = None
        neutral_prescribed_advection = None
        self.typeTi = IONS_T_I_NEGLECT
        if 'typeTi' in data:
            self.typeTi = int(data['typeTi'])
        if 'initial' in data:
            initial = data['initial']
        if 'prescribed' in data:
            prescribed = data['prescribed']
        if 'charged_prescribed_diffusion' in data:
            charged_prescribed_diffusion = data['charged_prescribed_diffusion']
        if 'neutral_prescribed_diffusion' in data:
            neutral_prescribed_diffusion = data['neutral_prescribed_diffusion']
        if 'charged_prescribed_advection' in data:
            charged_prescribed_advection = data['charged_prescribed_advection']
        if 'neutral_prescribed_advection' in data:
            neutral_prescribed_advection = data['neutral_prescribed_advection']
        if 'adv_interp_charged' in data:
            self.advectionInterpolationCharged.fromdict(data['adv_interp_charged'])
        if 'adv_interp_neutral' in data:
            self.advectionInterpolationNeutral.fromdict(data['adv_interp_neutral'])
        if 'initialTi' in data:
            initialTi = data['initialTi']
        iidx, pidx, spiidx, cpdidx, npdidx, cpaidx, npaidx = 0, 0, 0, 0, 0, 0, 0
        for i in range(len(Z)):
            if types[i] == IONS_PRESCRIBED:
                n = prescribed['x'][pidx:(pidx+Z[i]+1)]
                r = prescribed['r']
                t = prescribed['t']
                pidx += Z[i]+1
            else:
                n = initial['x'][iidx:(iidx+Z[i]+1)]
                r = initial['r']
                t = None #initial['t']
                iidx += Z[i]+1
            if self.typeTi==IONS_T_I_INCLUDE and initialTi is not None:
                T = initialTi['x'][i]
            else: 
                T = None
            if SPIMolarFraction[spiidx]>=0:
                SPIMolarFractionSingleSpecies = SPIMolarFraction[spiidx:spiidx+nShard]
                spiidx+=nShard
            else:
                SPIMolarFractionSingleSpecies = SPIMolarFraction[spiidx]
                spiidx+=1
            tritium = (names[i] in tritiumnames)
            hydrogen = (names[i] in hydrogennames)
            
            if charged_diffusion_modes[i] == ION_CHARGED_DIFFUSION_MODE_PRESCRIBED:
                cpd = charged_prescribed_diffusion['x'][cpdidx:(cpdidx+Z[i])]
                rcpd = charged_prescribed_diffusion['r']
                tcpd = charged_prescribed_diffusion['t']
                cpdidx += Z[i]
            else:
                cpd=None
                rcpd=None
                tcpd=None
                
            if neutral_diffusion_modes[i] == ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED:
                npd = neutral_prescribed_diffusion['x'][npdidx:(npdidx+1)]
                rnpd = neutral_prescribed_diffusion['r']
                tnpd = neutral_prescribed_diffusion['t']
                npdidx += 1
            else:
                npd=None
                rnpd=None
                tnpd=None
                
            if charged_advection_modes[i] == ION_CHARGED_ADVECTION_MODE_PRESCRIBED:
                cpa = charged_prescribed_advection['x'][cpaidx:(cpaidx+Z[i])]
                rcpa = charged_prescribed_advection['r']
                tcpa = charged_prescribed_advection['t']
                cpaidx += Z[i]
            else:
                cpa=None
                rcpa=None
                tcpa=None
                
            if neutral_advection_modes[i] == ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED:
                npa = neutral_prescribed_advection['x'][npaidx:(npaidx+1)]
                rnpa = neutral_prescribed_advection['r']
                tnpa = neutral_prescribed_advection['t']
                npaidx += 1
            else:
                npa=None
                rnpa=None
                tnpa=None

            self.addIon(name=names[i], Z=Z[i], isotope=isotopes[i], SPIMolarFraction=SPIMolarFractionSingleSpecies, iontype=types[i], opacity_mode=opacity_modes[i], 
                charged_diffusion_mode=charged_diffusion_modes[i], charged_prescribed_diffusion = cpd, rChargedPrescribedDiffusion=rcpd, tChargedPrescribedDiffusion = tcpd,
                neutral_diffusion_mode=neutral_diffusion_modes[i], neutral_prescribed_diffusion = npd, rNeutralPrescribedDiffusion=rnpd, tNeutralPrescribedDiffusion = tnpd,
                charged_advection_mode=charged_advection_modes[i], charged_prescribed_advection = cpa, rChargedPrescribedAdvection=rcpa, tChargedPrescribedAdvection = tcpa,
                neutral_advection_mode=neutral_advection_modes[i], neutral_prescribed_advection = npa, rNeutralPrescribedAdvection=rnpa, tNeutralPrescribedAdvection = tnpa,
                T=T, n=n, r=r, t=t, tritium=tritium, hydrogen=hydrogen)

        if 'ionization' in data:
            self.ionization = int(data['ionization'])

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this Ions object.
        """

        Z       = self.getCharges()
        itypes  = self.getTypes()
        iopacity_modes =self.getOpacityModes()
        icharged_diffusion_modes =self.getChargedDiffusionModes()
        ineutral_diffusion_modes =self.getNeutralDiffusionModes()
        icharged_advection_modes =self.getChargedAdvectionModes()
        ineutral_advection_modes =self.getNeutralAdvectionModes()
        isotopes     = self.getIsotopes()
        initial = None
        initialTi = None
        prescribed = None
        sourceterm = None
        sourceterm_types = []
        charged_prescribed_diffusion = None
        neutral_prescribed_diffusion = None
        charged_prescribed_advection = None
        neutral_prescribed_advection = None
        names   = ""
        init_equil = []
        initialNi = []

        hydrogennames = ""
        tritiumnames = ""

        SPIMolarFraction = None

        for ion in self.ions:
            names += '{};'.format(ion.getName())

            if ion.tritium:
                tritiumnames += '{};'.format(ion.getName())
            elif ion.hydrogen:
                hydrogennames += '{};'.format(ion.getName())

            # Set prescribed/initial density
            init_equil.append(1 if ion.initializeToEquilibrium() else 0)
            if ion.ttype != IONS_PRESCRIBED:
                ni = ion.getDensity()
                if ni is None:
                    ni = np.zeros((ion.Z+1, self.r.size))

                if initial is None:
                    initial = np.copy(ni)
                else:
                    initial = np.concatenate((initial, ni))
            else:
                if prescribed is None:
                    prescribed = np.copy(ion.getDensity())
                else:
                    prescribed = np.concatenate((prescribed, ion.getDensity()))

            # Construct source term
            sourceterm_types.append(ion.getSourceType())
            if sourceterm is None:
                sourceterm = np.copy(ion.getSourceDensity())
            else:
                n1 = sourceterm
                n2 = ion.getSourceDensity()

                if n1.shape[1] != n2.shape[1]:
                    if n1.shape[1] == 1:
                        n1 = repmat(n1[:,0].reshape((n1.shape[0],1)), 1, n2.shape[1])
                    elif n2.shape[1] == 1:
                        n2 = repmat(n2[:,0].reshape((n2.shape[0],1)), 1, n1.shape[1])
                    else:
                        raise EquationException("All ion sources must be defined in the same time points.")

                sourceterm = np.concatenate((n1, n2))

            if initialTi is None:
                initialTi = np.copy(ion.getTemperature())
            else:
                initialTi = np.concatenate((initialTi, ion.getTemperature()))

            initialNi.append(ion.getInitialSpeciesDensity())
                
            if SPIMolarFraction is None:
                SPIMolarFraction = np.copy(ion.getSPIMolarFraction())
            else:
                SPIMolarFraction = np.concatenate((SPIMolarFraction, ion.getSPIMolarFraction()))
                
            if ion.getChargedDiffusionMode()==ION_CHARGED_DIFFUSION_MODE_PRESCRIBED:
                if charged_prescribed_diffusion is None:
                    charged_prescribed_diffusion = np.copy(ion.getChargedPrescribedDiffusion())
                else:
                    charged_prescribed_diffusion = np.concatenate((charged_prescribed_diffusion, ion.getChargedPrescribedDiffusion()))
           
            if ion.getNeutralDiffusionMode()==ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED:
                if neutral_prescribed_diffusion is None:
                    neutral_prescribed_diffusion = np.copy(ion.getNeutralPrescribedDiffusion())
                else:
                    neutral_prescribed_diffusion = np.concatenate((neutral_prescribed_diffusion, ion.getNeutralPrescribedDiffusion()))

            if ion.getChargedAdvectionMode()==ION_CHARGED_ADVECTION_MODE_PRESCRIBED:
                if charged_prescribed_advection is None:
                    charged_prescribed_advection = np.copy(ion.getChargedPrescribedAdvection())
                else:
                    charged_prescribed_advection = np.concatenate((charged_prescribed_advection, ion.getChargedPrescribedAdvection()))
           
            if ion.getNeutralAdvectionMode()==ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED:
                if neutral_prescribed_advection is None:
                    neutral_prescribed_advection = np.copy(ion.getNeutralPrescribedAdvection())
                else:
                    neutral_prescribed_advection = np.concatenate((neutral_prescribed_advection, ion.getNeutralPrescribedAdvection()))
                
        data = {
            'names': names,
            'Z': Z,
            'isotopes':isotopes,
            'SPIMolarFraction':SPIMolarFraction,
            'types': itypes,
            'opacity_modes':iopacity_modes,
            'charged_diffusion_modes':icharged_diffusion_modes,
            'neutral_diffusion_modes':ineutral_diffusion_modes,
            'charged_advection_modes':icharged_advection_modes,
            'neutral_advection_modes':ineutral_advection_modes
        }

        if len(tritiumnames) > 0:
            data['tritiumnames'] = tritiumnames
        if len(hydrogennames) > 0:
            data['hydrogennames'] = hydrogennames

        if initial is not None and len(initial) > 0:
            for i in range(len(initial)):
                if initial[i] is None:
                    initial[i] = np.zeros((self.Z+1, self.r.size,))
            initial = np.array(initial)

            data['initial'] = {
                'r': self.r,
                'x': initial
            }

        if prescribed is not None:
            data['prescribed'] = {
                'r': self.r,
                't': self.t,
                'x': prescribed
            }

        if charged_prescribed_diffusion is not None:
            data['charged_prescribed_diffusion'] = {
                'r': self.rChargedPrescribedDiffusion,
                't': self.tChargedPrescribedDiffusion,
                'x': charged_prescribed_diffusion
            }
            
        if neutral_prescribed_diffusion is not None:
            data['neutral_prescribed_diffusion'] = {
                'r': self.rNeutralPrescribedDiffusion,
                't': self.tNeutralPrescribedDiffusion,
                'x': neutral_prescribed_diffusion
            }
            
        if charged_prescribed_advection is not None:
            data['charged_prescribed_advection'] = {
                'r': self.rChargedPrescribedAdvection,
                't': self.tChargedPrescribedAdvection,
                'x': charged_prescribed_advection
            }
            
        if neutral_prescribed_advection is not None:
            data['neutral_prescribed_advection'] = {
                'r': self.rNeutralPrescribedAdvection,
                't': self.tNeutralPrescribedAdvection,
                'x': neutral_prescribed_advection
            }

        if self.tSourceTerm is not None:
            data['ion_source_types'] = sourceterm_types
            data['ion_source'] = {
                't': self.tSourceTerm,
                'x': sourceterm
            }
        
        # Flux limiter settings
        data['adv_interp_charged'] = self.advectionInterpolationCharged.todict()
        data['adv_interp_neutral'] = self.advectionInterpolationNeutral.todict()
            
        data['initialTi'] = {
            'r': self.r,
            'x': initialTi
        }
        data['ionization'] = self.ionization
        data['typeTi'] = self.typeTi
        
        # Initial equilibrium
        for i in range(len(initialNi)):
            if initialNi[i] is None:
                initialNi[i] = np.zeros((self.r.size,))
        initialNi = np.array(initialNi)

        data['init_equilibrium'] = init_equil
        data['initialNi'] = {
            'r': self.r,
            'x': initialNi
        }

        return data
            

    def verifySettings(self):
        """
        Verify that all settings are consistent.
        """
        # Make sure there are no double names
        for i in range(0, len(self.ions)):
            for j in range(0, len(self.ions)):
                if i == j: continue

                if self.ions[i].getName() == self.ions[j].getName():
                    raise EquationException("ions: More than one ion species is named '{}'.".format(self.ions[i].getName()))
            
            self.ions[i].verifySettings()
        
        if (self.ionization != IONIZATION_MODE_FLUID) and (self.ionization != IONIZATION_MODE_KINETIC) and (self.ionization != IONIZATION_MODE_KINETIC_APPROX_JAC):
            raise EquationException("ions: Invalid ionization mode: {}.".format(self.ionization))
 

    def getFreeElectronDensity(self, t=0):
        """
        Returns the plasma free electron density at the given time index, based
        on the prescribed/initialized ion densities.

        :param int t: Index of time for which to retrieve the free electron density.
        """
        n_free = np.zeros( self.r.shape )

        for ion in self.ions:
            for Z0 in range(1,ion.Z + 1):
                if len( ion.n.shape ) == 3:
                    n_free = n_free + Z0 * ion.n[Z0,t,:]
                elif len( ion.n.shape ) == 2:
                    n_free = n_free + Z0 * ion.n[Z0,:]
                
        return n_free, self.r


