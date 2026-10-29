
import numpy as np
import scipy.interpolate

from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.IonSpecies import IonSpecies, IONS_PRESCRIBED, IONIZATION_MODE_FLUID, IONIZATION_MODE_KINETIC, IONIZATION_MODE_KINETIC_APPROX_JAC, IONIZATION_MODE_FLUID_RE, ION_OPACITY_MODE_TRANSPARENT, ION_CHARGED_DIFFUSION_MODE_NONE, ION_CHARGED_DIFFUSION_MODE_PRESCRIBED, ION_NEUTRAL_DIFFUSION_MODE_NONE, ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED, ION_CHARGED_ADVECTION_MODE_NONE, ION_CHARGED_ADVECTION_MODE_PRESCRIBED, ION_NEUTRAL_ADVECTION_MODE_NONE, ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED, ION_SOURCE_NONE, ION_SOURCE_PRESCRIBED, ION_SOURCE_PRESCRIBED_VOLUMETRIC
from . UnknownQuantity import UnknownQuantity
from .. import AdvectionInterpolation
from ... DREAMException import DREAMException
from ... helpers import scal
from .IonSpecies import IONS_DYNAMIC
from .. import MolecularReactionData as reaction_data


# Model to use for ion heat
IONS_T_I_NEGLECT = 1
IONS_T_I_INCLUDE = 2
IONS_T_N_NEGLECT = 1
IONS_T_N_INCLUDE = 2

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
        self.rSourceTerm                 = None

        self.ionization = ionization
        self.typeTi = IONS_T_I_NEGLECT
        self.reioniz_scale = 1

        self.advectionInterpolationCharged = AdvectionInterpolation.AdvectionInterpolation(kinetic=False)
        self.advectionInterpolationNeutral = AdvectionInterpolation.AdvectionInterpolation(kinetic=False)

        self.reactions_enabled = False
        self.enabled_reaction_names = []
        self.user_rates = {}
        self.typeTn = IONS_T_N_NEGLECT

    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, Z0=None, isotope=0, SPIMolarFraction=-1, opacity_mode=ION_OPACITY_MODE_TRANSPARENT,
        charged_diffusion_mode=ION_CHARGED_DIFFUSION_MODE_NONE, charged_prescribed_diffusion=None, rChargedPrescribedDiffusion=None, tChargedPrescribedDiffusion=None,
        neutral_diffusion_mode=ION_NEUTRAL_DIFFUSION_MODE_NONE, neutral_prescribed_diffusion=None, rNeutralPrescribedDiffusion=None, tNeutralPrescribedDiffusion=None,
        charged_advection_mode=ION_CHARGED_ADVECTION_MODE_NONE, charged_prescribed_advection=None, rChargedPrescribedAdvection=None, tChargedPrescribedAdvection=None,
        neutral_advection_mode=ION_NEUTRAL_ADVECTION_MODE_NONE, neutral_prescribed_advection=None, rNeutralPrescribedAdvection=None, tNeutralPrescribedAdvection=None,
        t_transp_expdecay_all_cs = None, t_transp_start_expdecay_all_cs = 0, diffusion_initial_all_cs = None, diffusion_final_all_cs = 0, diffusion_offset_all_cs = 0, advection_initial_all_cs = None, advection_final_all_cs = 0, advection_offset_all_cs = 0, 
        r_expdecay_all_cs = None, t_expdecay_all_cs = None, init_equil=False, T=None, n=None, r=None, t=None, tritium=False, hydrogen=False, Tn = None):

        """
        Adds a new ion species to the plasma.

        :param str name:        Name by which the ion species will be referred to.
        :param int Z:           Ion charge number (and summed nuclear charges for molecular species).
        :param int isotope:            Ion mass number.
        :param int iontype:     Method to use for evolving ions in time.
        :param int Z0:          Charge state to populate (used for populating exactly one charge state for the ion).
        :param n:               Ion density (can be either a scalar, 1D array or 2D array, depending on the other input parameters)
        :param float SPIMolarFraction: Molar fraction of the SPI injection (if any). A negative value means that this species is not part of the SPI injection
        :param numpy.ndarray r: Radial grid on which the input density is defined.
        :param T:               Ion initial temperature (can be scalar for uniform temperature, otherwise 1D array matching `r` in size)
        :param Tn:              Ion initial neutral temperature (can be scalar for uniform temperature, otherwise 1D array matching `r` in size)
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
            diffusion_offset_all_cs = diffusion_offset_all_cs, advection_initial_all_cs = advection_initial_all_cs, 
            advection_final_all_cs = advection_final_all_cs, advection_offset_all_cs = advection_offset_all_cs, 
            r_expdecay_all_cs = r_expdecay_all_cs, t_expdecay_all_cs = t_expdecay_all_cs,
            init_equil=init_equil, T=T, n=n, r=r, t=t, interpr=self.r, interpt=None, tritium=tritium, hydrogen=hydrogen, Tn=Tn)

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


    def addIonSource(self, species, dNdt=None, t=None, Z0=0, source_type=ION_SOURCE_PRESCRIBED, r =None):
        """
        Add a source term for the specified ion species.

        :param species: Name of the ion species to add this source term to.
        :param dNdt:    Number of particles to add per unit time. Can be:
                - Scalar: constant in space and time for charge state Z0
                - 1D array with t: time-dependent for charge state Z0
                - 1D array with r: radially-dependent for charge state Z0
                - 2D array (Z+1, nt) with t: time-dependent for ALL charge states (boundary source)
                - 3D array (Z+1, nt, nr) with t and r: time- and radially-dependent for ALL charge states (volumetric source)
        :param t:       Time grid associated with ``dNdt`` (if any).
        :param r:       Radial grid associated with ``dNdt`` (if any).
        :param Z0:      For scalar or 1D ``dNdt``, the charge state for which to add the source term.
        :param source_type: Type of source term to add (e.g., prescribed, constant).
        """
        if t is None:
            t = self.tSourceTerm
        elif self.tSourceTerm is not None and not np.all(t == self.tSourceTerm):
            raise EquationException("The time grid used for ion sources must be the same for all ion species.")
        
        if r is None:
            r = self.rSourceTerm

        found = False
        for ion in self.ions:
            if ion.name == species:
                ion.initialize_source(n=dNdt, t=t, Z0=Z0, source_type=source_type, r=r)
                found = True

                self.tSourceTerm = ion.getSourceTime()
                if ion.getSourceRadialGrid() is not None:
                    self.rSourceTerm = ion.getSourceRadialGrid()

        if not found:
            raise EquationException(f"No ion species with name '{species}' has been added to the simulation. Unable to add source term.")

    def enableReactions(self, selection=None):
        """
        Enable molecular reactions for the ion species.
        This will call a method to add all species that are required for the
        selected reactions to the list of ions. If selection is None,
        all (possible) reactions will be enabled.
        Otherwise, only the selected reactions will be enabled."""

        reactions = list(reaction_data.MOLECULAR_REACTIONS)

        #If we select certain reactions, filter the list of reactions to only include those.
        if selection is not None:
            selected = []
            for process, requests in selection.items():
                candidates = [
                                r for r in reactions
                                if r["classification"] == process
                            ]

                if not candidates:
                    raise EquationException(
                        f"Unknown or empty reaction process '{process}'."
                    )

                if requests == "all":
                    selected.extend(candidates)
                    continue


                for request in requests:
                    # A list preserves repeated reactants.
                    wanted = sorted(tuple(x) for x in request["reactants"])
                    matches = []

                    for reaction in candidates:
                        actual = []

                        for name, Z0, coefficient in reaction["reactants"]:
                            actual.extend([(name, Z0)] * coefficient)

                        if sorted(actual) == wanted:
                            matches.append(reaction)

                    if not matches:
                        raise EquationException(
                            f"No reaction matches {request}."
                        )

                    for reaction in matches:
                        if reaction not in selected:
                            selected.append(reaction)

            reactions = selected

        self.addRequiredSpecies(reactions)

        self.reactions_enabled = True
        self.enabled_reaction_names = []
        species_names = []
        for ion in self.ions:
            species_names.append(ion.getName())

        # Check each reaction.
        for reaction in reactions:
            all_reactants_exist = True

            for reactant in reaction["reactants"]:
                name = reactant[0]

                # Electrons do not need to be declared as ion species.
                if name in ("e", "RE"):
                    continue

                if name not in species_names:
                    all_reactants_exist = False
                    print(f"Warning: Reaction '{reaction['name']}' is disabled because reactant species '{name}' does not exist.")
                    break

            # Save the reaction's name if all its reactant species exist.
            if all_reactants_exist:
                self.enabled_reaction_names.append(reaction["name"])

    def addRequiredSpecies(self, reactions):
        """
        Add all species required by reactions with available reactants.
        """
        charges = {}
        initial_names = set()

        for ion in self.ions:
            name = ion.getName()
            Z = ion.getZ()

            charges[name] = Z
            initial_names.add(name)

        changed = True

        while changed:
            changed = False ##Assume we will not discove somet

            for reaction in reactions:
                all_reactants_exist = True

                for name, Z0, coefficient in reaction["reactants"]:
                    if name in ("e", "RE"):
                        continue

                    if name not in charges: #If the name is not in charges we add it in the next loop, and then we have to go though it again
                        all_reactants_exist = False
                        break

                if not all_reactants_exist:
                    continue

                for name, Z0, coefficient in reaction["products"]:
                    if name in ("e", "RE"):
                        continue

                    if name in charges:
                        continue

                    if name not in reaction_data.SPECIES_Z:
                        raise EquationException(
                            f"Add species metadata for '{name}' in data.py."
                        )

                    charges[name] = reaction_data.SPECIES_Z[name]
                    changed = True

        # Check charge states before adding species to self.ions.
        for reaction in reactions:
            all_reactants_exist = True

            for name, Z0, coefficient in reaction["reactants"]:
                if name in ("e", "RE"):
                    continue

                if name not in charges:
                    all_reactants_exist = False
                    break

            if not all_reactants_exist:
                continue

            participants = reaction["reactants"] + reaction["products"]

            for name, Z0, coefficient in participants:
                if name in ("e", "RE"):
                    continue

                if Z0 < 0 or Z0 > charges[name]:
                    raise EquationException(
                        f"Unsupported charge state {Z0} for '{name}'."
                    )

        # Add the newly discovered species with zero density.
        for name in charges:
            if name not in initial_names:
                self.addIon(
                    name=name,
                    Z=charges[name],
                    iontype=IONS_DYNAMIC,
                    Z0=0,
                    n=0,
                    r=self.r,
                )

    def setReactionRate(self, name, rate, T = None, n=None):
        """
        Replace the pre-calculated rate of a molecular reaction with
        user-supplied data.

        :param str name:  Name of the reaction (as in MolecularReactionData.py).
        :param T:         Temperature grid [eV], strictly increasing.
                          If not given, the rate is independent of temperature.
        :param rate:      Rate coefficient <sigma v> [m^3/s], shape (nT,) if
                          n is None, otherwise shape (nT, nn).
        :param n:         Density grid [m^-3], strictly increasing.
                          If not given, the rate is independent of density.
        """
        names = [r["name"] for r in reaction_data.MOLECULAR_REACTIONS]
        if name not in names:
            raise EquationException(f"Unknown molecular reaction '{name}'.")


        rate = np.asarray(rate, dtype=float)

        if T is None:
            T = np.array([1e-2, 1e5])
            rate = np.stack((rate, rate))
        if n is None:
            n = np.array([1e14, 1e22])
            rate = np.stack((rate, rate), axis=-1)
        T = np.asarray(T, dtype=float)
        n = np.asarray(n, dtype=float)

        if rate.shape != (T.size, n.size):
            raise EquationException(f"Rate '{name}': 'rate' must have shape (nT, nn) = ({T.size}, {n.size}), but has shape {rate.shape}.")
        if T.size < 2 or np.any(np.diff(T) <= 0) or np.any(np.diff(n) <= 0):
            raise EquationException(f"Rate '{name}': 'T' and 'n' must be strictly increasing with at least two points.")
        if np.any(T <= 0) or np.any(n <= 0) or np.any(rate <= 0):
            raise EquationException(f"Rate '{name}': 'T', 'n' and 'rate' must be positive (interpolation is done in log10).")

        self.user_rates[name] = {'n': n, 'T': T, 'rate': rate}




    def changeRadialGrid(self, r):
        """
        Change the radial grid used for the ion species.
        """
        for ion in self.ions:
            if ion.r.size == 1:
                ion.n = ion.n * np.ones(ion.n.shape[:-1] + (r.size,))
                ion.T = ion.T * np.ones(ion.T.shape[:-1] + (r.size,))
                ion.r = r
                ion.Tn = ion.Tn * np.ones(ion.Tn.shape[:-1] + (r.size,))
            else:
                fn = scipy.interpolate.interp1d(ion.r, ion.n, axis=-1, bounds_error=False, fill_value='extrapolate')
                fT = scipy.interpolate.interp1d(ion.r, ion.T, axis=-1, bounds_error=False, fill_value='extrapolate')
                fTn = scipy.interpolate.interp1d(ion.r, ion.Tn, axis=-1, bounds_error=False, fill_value='extrapolate')
                ion.n = fn(r)
                ion.T = fT(r)
                ion.Tn = fTn(r)
                ion.r = r

        self.r = r

    def setNeutralTemperatureEnabled(self, enabled=True):
        self.typeTn = (IONS_T_N_INCLUDE if enabled else IONS_T_N_NEGLECT)


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


    def setIonization(self, ionization=IONIZATION_MODE_FLUID, reioniz_scale=1.0):
        """
        Sets which model to use for ionization.

        :param int ionization: Flag indicating which model to use for ionization.
        """
        self.ionization=ionization

        if ionization == IONIZATION_MODE_FLUID_RE:
            self.reioniz_scale = reioniz_scale


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


    def setChargedDiffusion(
        self, species, mode, Drr=None, r=None, t=None,
        t_transp_expdecay_all_cs=None, t_transp_start_expdecay_all_cs=0,
        Drr_0=None, Drr_f=0,
        r_expdecay_all_cs=None, t_expdecay_all_cs=None,
        interpr=None, interpt=None
    ):
        """
        Set ion radial diffusion for charged particles.

        :param species:                        Species to apply transport to.
        :param mode:                           Type of transport to prescribe.
        :param Drr:                            Diffusion coefficient to prescribe.
        :param r:                              Radial grid on which ``Drr`` is given.
        :param t:                              Time grid on which ``Drr`` is given.
        :param t_transp_expdecay_all_cs:       e-folding time of transport coefficient decay (exponential decay).
        :param t_transp_start_expdecay_all_cs: Start time of exponential decay.
        :param Drr_0:                          Initial value of diffusion coefficient when decaying exponentially.
        :param Drr_f:                          Final value of diffusion coefficient when decaying exponentially.
        :param r_expdecay_all_cs:              Radial grid on which the coefficient should be defined.
        :param t_expdecay_all_cs:              Time grid on which the coefficient should be defined.
        :param interpr:                        Radial grid onto which ion transport coefficients should be interpolated.
        :param interpt:                        Time grid onto which ion transport coefficients should be interpolated.
        """
        ion = self.getIon(self.getIndex(species))
        ion.setChargedDiffusion(
            mode=mode, Drr=Drr, r=r, t=t,
            t_transp_expdecay_all_cs=t_transp_expdecay_all_cs,
            t_transp_start_expdecay_all_cs=t_transp_start_expdecay_all_cs,
            Drr_0=Drr_0, Drr_f=Drr_f, r_expdecay_all_cs=r_expdecay_all_cs,
            t_expdecay_all_cs=t_expdecay_all_cs,
            interpr=interpr, interpt=interpt
        )

        if mode == ION_CHARGED_DIFFUSION_MODE_PRESCRIBED:
            self.rChargedPrescribedDiffusion = ion.getRChargedPrescribedDiffusion()
            self.tChargedPrescribedDiffusion = ion.getTChargedPrescribedDiffusion()


    def setNeutralDiffusion(
        self, species, mode, Drr=None, r=None, t=None,
        t_transp_expdecay_all_cs=None, t_transp_start_expdecay_all_cs=0,
        Drr_0=None, Drr_f=0,
        r_expdecay_all_cs=None, t_expdecay_all_cs=None,
        interpr=None, interpt=None
    ):
        """
        Set ion radial diffusion for neutral particles.

        :param species:                        Species to apply transport to.
        :param mode:                           Type of transport to prescribe.
        :param Drr:                            Diffusion coefficient to prescribe.
        :param r:                              Radial grid on which ``Drr`` is given.
        :param t:                              Time grid on which ``Drr`` is given.
        :param t_transp_expdecay_all_cs:       e-folding time of transport coefficient decay (exponential decay).
        :param t_transp_start_expdecay_all_cs: Start time of exponential decay.
        :param Drr_0:                          Initial value of diffusion coefficient when decaying exponentially.
        :param Drr_f:                          Final value of diffusion coefficient when decaying exponentially.
        :param r_expdecay_all_cs:              Radial grid on which the coefficient should be defined.
        :param t_expdecay_all_cs:              Time grid on which the coefficient should be defined.
        :param interpr:                        Radial grid onto which ion transport coefficients should be interpolated.
        :param interpt:                        Time grid onto which ion transport coefficients should be interpolated.
        """
        ion = self.getIon(self.getIndex(species))
        ion.setNeutralDiffusion(
            mode=mode, Drr=Drr, r=r, t=t,
            t_transp_expdecay_all_cs=t_transp_expdecay_all_cs,
            t_transp_start_expdecay_all_cs=t_transp_start_expdecay_all_cs,
            Drr_0=Drr_0, Drr_f=Drr_f, r_expdecay_all_cs=r_expdecay_all_cs,
            t_expdecay_all_cs=t_expdecay_all_cs,
            interpr=interpr, interpt=interpt
        )

        if mode == ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED:
            self.rNeutralPrescribedDiffusion = ion.getRNeutralPrescribedDiffusion()
            self.tNeutralPrescribedDiffusion = ion.getTNeutralPrescribedDiffusion()


    def setChargedAdvection(
        self, species, mode, Ar=None, r=None, t=None,
        t_transp_expdecay_all_cs=None, t_transp_start_expdecay_all_cs=0,
        Ar_0=None, Ar_f=0,
        r_expdecay_all_cs=None, t_expdecay_all_cs=None,
        interpr=None, interpt=None
    ):
        """
        Set ion radial advection for charged particles.

        :param species:                        Species to apply transport to.
        :param mode:                           Type of transport to prescribe.
        :param Ar:                             Advection coefficient to prescribe.
        :param r:                              Radial grid on which ``Ar`` is given.
        :param t:                              Time grid on which ``Ar`` is given.
        :param t_transp_expdecay_all_cs:       e-folding time of transport coefficient decay (exponential decay).
        :param t_transp_start_expdecay_all_cs: Start time of exponential decay.
        :param Ar_0:                           Initial value of advection coefficient when decaying exponentially.
        :param Ar_f:                           Final value of advection coefficient when decaying exponentially.
        :param r_expdecay_all_cs:              Radial grid on which the coefficient should be defined.
        :param t_expdecay_all_cs:              Time grid on which the coefficient should be defined.
        :param interpr:                        Radial grid onto which ion transport coefficients should be interpolated.
        :param interpt:                        Time grid onto which ion transport coefficients should be interpolated.
        """
        ion = self.getIon(self.getIndex(species))
        ion.setChargedAdvection(
            mode=mode, Ar=Ar, r=r, t=t,
            t_transp_expdecay_all_cs=t_transp_expdecay_all_cs,
            t_transp_start_expdecay_all_cs=t_transp_start_expdecay_all_cs,
            Ar_0=Ar_0, Ar_f=Ar_f, r_expdecay_all_cs=r_expdecay_all_cs,
            t_expdecay_all_cs=t_expdecay_all_cs,
            interpr=interpr, interpt=interpt
        )

        if mode == ION_CHARGED_ADVECTION_MODE_PRESCRIBED:
            self.rChargedPrescribedAdvection = ion.getRChargedPrescribedAdvection()
            self.tChargedPrescribedAdvection = ion.getTChargedPrescribedAdvection()


    def setNeutralAdvection(
        self, species, mode, Ar=None, r=None, t=None,
        t_transp_expdecay_all_cs=None, t_transp_start_expdecay_all_cs=0,
        Ar_0=None, Ar_f=0,
        r_expdecay_all_cs=None, t_expdecay_all_cs=None,
        interpr=None, interpt=None
    ):
        """
        Set ion radial advection for neutral particles.

        :param species:                        Species to apply transport to.
        :param mode:                           Type of transport to prescribe.
        :param Ar:                             Advection coefficient to prescribe.
        :param r:                              Radial grid on which ``Ar`` is given.
        :param t:                              Time grid on which ``Ar`` is given.
        :param t_transp_expdecay_all_cs:       e-folding time of transport coefficient decay (exponential decay).
        :param t_transp_start_expdecay_all_cs: Start time of exponential decay.
        :param Ar_0:                           Initial value of advection coefficient when decaying exponentially.
        :param Ar_f:                           Final value of advection coefficient when decaying exponentially.
        :param r_expdecay_all_cs:              Radial grid on which the coefficient should be defined.
        :param t_expdecay_all_cs:              Time grid on which the coefficient should be defined.
        :param interpr:                        Radial grid onto which ion transport coefficients should be interpolated.
        :param interpt:                        Time grid onto which ion transport coefficients should be interpolated.
        """
        ion = self.getIon(self.getIndex(species))
        ion.setNeutralAdvection(
            mode=mode, Ar=Ar, r=r, t=t,
            t_transp_expdecay_all_cs=t_transp_expdecay_all_cs,
            t_transp_start_expdecay_all_cs=t_transp_start_expdecay_all_cs,
            Ar_0=Ar_0, Ar_f=Ar_f, r_expdecay_all_cs=r_expdecay_all_cs,
            t_expdecay_all_cs=t_expdecay_all_cs,
            interpr=interpr, interpt=interpt
        )

        if mode == ION_CHARGED_ADVECTION_MODE_PRESCRIBED:
            self.rNeutralPrescribedAdvection = ion.getRNeutralPrescribedAdvection()
            self.tNeutralPrescribedAdvection = ion.getTNeutralPrescribedAdvection()
    

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
        initialTn  = None
        charged_prescribed_diffusion = None
        neutral_prescribed_diffusion = None
        charged_prescribed_advection = None
        neutral_prescribed_advection = None
        ion_source_t = None
        ion_source_x = None
        ion_source_r = None
        self.typeTi = IONS_T_I_NEGLECT
        self.typeTn = IONS_T_N_NEGLECT
        if 'typeTi' in data:
            self.typeTi = int(scal(data['typeTi']))
        if 'typeTn' in data:
            self.typeTn = int(scal(data['typeTn']))
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
        if 'initialTn' in data:
            initialTn = data['initialTn']
        if 'init_equilibrium' in data:
            init_equilibrium = data['init_equilibrium']
        if 'initialNi' in data:
            initialNi = data['initialNi']

        if 'ion_source_types' in data:
            ion_source_types = data['ion_source_types']
            if len(ion_source_types) == 0:
                ion_source_types = None
        else:
            ion_source_types = None

        if 'ion_source' in data:
            ion_source_t = data['ion_source']['t']
            ion_source_x = data['ion_source']['x']
            
            if len(ion_source_x.shape) == 3:
                ion_source_r = data['ion_source'].get('r', None)
                if ion_source_r is None:
                    raise EquationException("The legacy 3D ion source ('ion_source') requires an explicit radial grid 'r'.")

        if 'ion_source_volumetric' in data:
            ion_source_t = data['ion_source_volumetric']['t']
            ion_source_x = data['ion_source_volumetric']['x']
            ion_source_r = data['ion_source_volumetric'].get('r', None)
            if ion_source_r is None:
                raise EquationException("The volumetric ion source ('ion_source_volumetric') requires an explicit radial grid 'r'.")


        iidx, pidx, spiidx, cpdidx, npdidx, cpaidx, npaidx, srcidx = 0, 0, 0, 0, 0, 0, 0, 0
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
            if self.typeTn==IONS_T_N_INCLUDE and initialTn is not None:
                Tn = initialTn['x'][i]
            else:
                Tn = None
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

            init_equil = (init_equilibrium[i] != 0)
            dens = initialNi['x'][i] if init_equil else n

            self.addIon(
                name=names[i], Z=Z[i], isotope=isotopes[i], SPIMolarFraction=SPIMolarFractionSingleSpecies,
                iontype=types[i], opacity_mode=opacity_modes[i], 
                # Charged diffusion
                charged_diffusion_mode=charged_diffusion_modes[i], charged_prescribed_diffusion = cpd,
                rChargedPrescribedDiffusion=rcpd, tChargedPrescribedDiffusion = tcpd,
                # Neutral diffusion
                neutral_diffusion_mode=neutral_diffusion_modes[i], neutral_prescribed_diffusion = npd,
                rNeutralPrescribedDiffusion=rnpd, tNeutralPrescribedDiffusion = tnpd,
                # Charged advection
                charged_advection_mode=charged_advection_modes[i], charged_prescribed_advection = cpa,
                rChargedPrescribedAdvection=rcpa, tChargedPrescribedAdvection = tcpa,
                # Neutral advection
                neutral_advection_mode=neutral_advection_modes[i], neutral_prescribed_advection = npa,
                rNeutralPrescribedAdvection=rnpa, tNeutralPrescribedAdvection = tnpa,
                T=T, Tn = Tn,n=dens, r=r, t=t, tritium=tritium, hydrogen=hydrogen, init_equil=init_equil)

            # Load ion source
            if ion_source_types is not None and ion_source_types[i] == ION_SOURCE_PRESCRIBED:
                self.addIonSource(names[i], dNdt=ion_source_x[srcidx,:], t=ion_source_t, source_type=ION_SOURCE_PRESCRIBED)
                srcidx += Z[i] + 1  # Move to next ion's position
            elif ion_source_types is not None and ion_source_types[i] == ION_SOURCE_PRESCRIBED_VOLUMETRIC:
                
                # Extract the full 3D source data for this ion species
                if len(ion_source_x.shape) == 3:
                    dNdt_data = ion_source_x[srcidx:srcidx+(Z[i]+1), :, :]
                    srcidx += Z[i] + 1  # Move to next ion's position
                else:
                    raise EquationException(f"Unexpected ion_source shape for volumetric source: {ion_source_x.shape}")
                
                self.addIonSource(names[i], dNdt=dNdt_data, t=ion_source_t, r=ion_source_r, source_type=ION_SOURCE_PRESCRIBED_VOLUMETRIC)

        if 'ionization' in data:
            self.ionization = int(scal(data['ionization']))

        if 'reioniz_scale' in data:
            self.reioniz_scale = float(scal(data['reioniz_scale']))

        if "reactions" in data:
            reactions = data["reactions"]

            self.reactions_enabled = bool(scal(reactions.get("enabled", 0)))

            names = reactions.get("names", "")
            self.enabled_reaction_names = [
                name for name in names.split(";") if name
            ]

            if "user_rates" in reactions:
                ur = reactions["user_rates"]
                on, oT, ok = 0, 0, 0
                for i, name in enumerate(ur["names"].split(";")[:-1]):
                    nn, nT = int(ur["nn"][i]), int(ur["nT"][i])
                    self.user_rates[name] = {
                        'n': ur["n"][on:on+nn],
                        'T': ur["T"][oT:oT+nT],
                        'rate': ur["rate"][ok:ok+nT*nn].reshape((nT, nn))
                    }
                    on += nn
                    oT += nT
                    ok += nT*nn

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
        initialTn = None
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

                # Ensure both arrays have same number of dimensions
                if len(n1.shape) == 2 and len(n2.shape) == 3:
                    # Expand n1 from (Z+1, nt) to (Z+1, nt, 1)
                    n1 = n1[:, :, np.newaxis]
                elif len(n1.shape) == 3 and len(n2.shape) == 2:
                    # Expand n2 from (Z+1, nt) to (Z+1, nt, 1)
                    n2 = n2[:, :, np.newaxis]

                # Now match time dimensions
                if n1.shape[1] != n2.shape[1]:
                    if n1.shape[1] == 1:
                        n1 = np.repeat(n1, n2.shape[1], axis=1)
                    elif n2.shape[1] == 1:
                        n2 = np.repeat(n2, n1.shape[1], axis=1)
                    else:
                        raise EquationException("All ion sources must be defined in the same time points.")
                
                # Match radial dimensions (only if both are 3D)
                if len(n1.shape) == 3 and len(n2.shape) == 3:
                    if n1.shape[2] != n2.shape[2]:
                        if n1.shape[2] == 1:
                            n1 = np.repeat(n1, n2.shape[2], axis=2)
                        elif n2.shape[2] == 1:
                            n2 = np.repeat(n2, n1.shape[2], axis=2)
                        else:
                            raise EquationException("All ion sources must be defined on the same radial grid.")

                sourceterm = np.concatenate((n1, n2))

            if initialTi is None:
                initialTi = np.copy(ion.getTemperature())
            else:
                initialTi = np.concatenate((initialTi, ion.getTemperature()))
            if initialTn is None:
                initialTn = np.copy(ion.getNeutralTemperature())
            else:
                initialTn = np.concatenate(
                    (initialTn, ion.getNeutralTemperature()))

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

            if len(sourceterm.shape) == 2:
                data['ion_source'] = {
                    't': self.tSourceTerm,
                    'x': sourceterm
                }
            elif len(sourceterm.shape) == 3:
                data['ion_source_volumetric'] = {
                    't': self.tSourceTerm,
                    'r': self.rSourceTerm,
                    'x': sourceterm
                }
            else:
                raise EquationException(f"Unexpected ion source shape: {sourceterm.shape}.")


        # Flux limiter settings
        data['adv_interp_charged'] = self.advectionInterpolationCharged.todict()
        data['adv_interp_neutral'] = self.advectionInterpolationNeutral.todict()

        data['initialTi'] = {
            'r': self.r,
            'x': initialTi
        }
        data['initialTn'] = {
            'r': self.r,
            'x': initialTn
        }

        data['ionization'] = self.ionization
        data['typeTi'] = self.typeTi
        data['typeTn'] = self.typeTn
        data['reioniz_scale'] = self.reioniz_scale

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

        data["reactions"] = {
        "enabled": int(self.reactions_enabled),
        "names": ";".join(self.enabled_reaction_names),
            }

        if self.user_rates:
            data["reactions"]["user_rates"] = {
                "names": "".join(f"{name};" for name in self.user_rates),
                "nn": [r['n'].size for r in self.user_rates.values()],
                "nT": [r['T'].size for r in self.user_rates.values()],
                "n": np.concatenate([r['n'] for r in self.user_rates.values()]),
                "T": np.concatenate([r['T'] for r in self.user_rates.values()]),
                "rate": np.concatenate([r['rate'].flatten() for r in self.user_rates.values()]),
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

        if (self.ionization != IONIZATION_MODE_FLUID) and (self.ionization != IONIZATION_MODE_KINETIC) and (self.ionization != IONIZATION_MODE_KINETIC_APPROX_JAC) and (self.ionization != IONIZATION_MODE_FLUID_RE):
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
