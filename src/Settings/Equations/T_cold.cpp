/**
 * Definition of equations relating to the cold electron temperature.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/Fluid/CollisionalEnergyTransferKineticTerm.hpp"
#include "DREAM/Equations/Fluid/IonisationHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Equations/Fluid/OhmicHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "DREAM/Equations/Fluid/SPIHeatAbsorbtionTerm.hpp"
#include "DREAM/Equations/Fluid/IonSPIIonizLossTerm.hpp"
#include "DREAM/Equations/Fluid/ElectronHeatTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Equations/Fluid/HaloRegionHeatLossTerm.hpp"
#include "DREAM/Equations/Fluid/NBIElectronTerm.hpp"
#include "DREAM/NBIHandler.hpp"



using namespace DREAM;
using namespace std;


#define MODULENAME "eqsys/T_cold"
#define MODULENAME_SPI "eqsys/spi"
#define MODULENAME_ION "eqsys/n_i"
#define MODULENAME_NRE "eqsys/n_re"


/**
 * Define options for the electron temperature module.
 */
void SimulationGenerator::DefineOptions_T_cold(Settings *s) {
	DefineOptions_T_cold_inner(s, MODULENAME);

	// Possibility for switching equation
	DefineOptions_TriggerCondition(s, MODULENAME "/switch");

	// Also define settings for the 'alternative' equation
	DefineOptions_T_cold_inner(s, MODULENAME "/switch/equation");
}

void SimulationGenerator::DefineOptions_T_cold_inner(Settings *s, const string& modulename) {
    s->DefineSetting(modulename + "/type", "Type of equation to use for determining the electron temperature evolution", (int_t)OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED);
    s->DefineSetting(modulename + "/recombination", "Whether to include recombination radiation (true) or ionization energy loss (false)", (bool)false);
    s->DefineSetting(modulename + "/halo_region_losses", "Whether to include losses through the halo region (true) or not (false)", (bool)false);
    s->DefineSetting(modulename + "/include_NBI", "Whether to include NBI heating term in T_cold evolution", (bool)false);
    DefineOptions_T_cold_NBI(s, modulename);

    // Prescribed data (in radius+time)
    DefineDataRT(modulename, s, "data");

    // Prescribed initial profile (when evolving T self-consistently)
    DefineDataR(modulename, s, "init");

    // Transport settings
    DefineOptions_Transport(modulename, s, false);
}
/**
 * Define options for the NBI heating term.
 */
void SimulationGenerator::DefineOptions_T_cold_NBI(Settings *s, const string& modulename) {
    const string NBI_PATH = string(modulename) + "/NBI";

    s->DefineSetting(modulename + "/NBI/enabled", "Enable/disable NBI heating", (bool)false);
    // Beam geometry settings
    s->DefineSetting(modulename + "/NBI/s_max", "Max beamline length to integrate", (real_t)2);
    s->DefineSetting(modulename + "/NBI/r_beam", "Beam radius", (real_t)0.1);
    s->DefineSetting(modulename + "/NBI/P0", "Beam starting point (x,y,z)", 3, (real_t*)nullptr);
    s->DefineSetting(modulename + "/NBI/n",  "Beam direction vector", 3, (real_t*)nullptr);
    s->DefineSetting(modulename + "/NBI/energy_fractions", "Energy fractions for multi-energy components", 3, (real_t*)nullptr);

    // Beam physics settings
    s->DefineSetting(modulename + "/NBI/Ti_beam", "Thermal ion temperature [eV]", (real_t)4.8e-15);
    s->DefineSetting(modulename + "/NBI/m_i_beam", "Ion mass of beam [kg]", (real_t)Constants::mD);

    // Plasma and beam parameters
    s->DefineSetting(modulename + "/NBI/R0", "Major Radius", (real_t)1.0);
    s->DefineSetting(modulename + "/NBI/j_B/t", "Radi grid for beam current density", 0, (real_t *)nullptr);
    s->DefineSetting(modulename + "/NBI/j_B/x", "Beam current density values", 0, (real_t *)nullptr);
    s->DefineSetting(modulename + "/NBI/j_B/tinterp", "Interpolation method for j_B", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
    s->DefineSetting(modulename + "/NBI/gaussian_profile", "Gaussian profile type for NBI: 0=disabled, 1=TCV, 2=ITER", (int_t)0);
    s->DefineSetting(modulename + "/NBI/P_NBI/t", "Time dependant power", 0, (real_t *)nullptr);
    s->DefineSetting(modulename + "/NBI/P_NBI/x", "Power deposition profile", 0, (real_t *)nullptr);
    s->DefineSetting(modulename + "/NBI/P_NBI/tinterp", "Interpolation method for P_NBI", (int_t)OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
}


/**
 * Construct the equation for the cold electron temperature.
 */
void SimulationGenerator::ConstructEquation_T_cold(
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
	const len_t id_T_cold = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);

	// Construct main equation
	ConstructEquation_T_cold_inner(MODULENAME, id_T_cold, eqsys, s, adas, nist, amjuel, &oqty_terms->T_cold);

	enum OptionConstants::eqn_trigger_type switchtype =
		(enum OptionConstants::eqn_trigger_type)s->GetInteger(MODULENAME "/switch/condition");
	
	// Set alternative equation?
	if (switchtype != OptionConstants::EQN_TRIGGER_TYPE_NONE) {
		eqsys->SetAssignToAlternativeEquation(id_T_cold, true);

		ConstructEquation_T_cold_inner(MODULENAME "/switch/equation", id_T_cold, eqsys, s, adas, nist, amjuel, &oqty_terms->T_cold);

		EquationTriggerCondition *trig = LoadTriggerCondition(s, MODULENAME "/switch", eqsys->GetFluidGrid(), eqsys->GetUnknownHandler());
		eqsys->SetTriggerCondition(id_T_cold, trig);

		eqsys->SetAssignToAlternativeEquation(id_T_cold, false);
	}

	ConstructEquation_W_cold(eqsys, s);
}

void SimulationGenerator::ConstructEquation_T_cold_inner(
	const string& modulename, const len_t id_T_cold,
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
    struct OtherQuantityHandler::T_terms *oqty_terms
) {
    enum OptionConstants::uqty_T_cold_eqn type =
		(enum OptionConstants::uqty_T_cold_eqn)s->GetInteger(modulename + "/type");

	FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    len_t id_W_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_n_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    len_t id_j_ohm  = unknowns->GetUnknownID(OptionConstants::UQTY_J_OHM);

    switch (type) {
        case OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED:
            ConstructEquation_T_cold_prescribed(modulename, eqsys, s, id_T_cold);
            break;

        case OptionConstants::UQTY_T_COLD_SELF_CONSISTENT: {
            ConstructEquation_T_cold_selfconsistent(
				modulename, eqsys, s, adas, nist, amjuel, oqty_terms,
				id_T_cold, id_T_cold,
				id_W_cold, id_n_cold, id_j_ohm
			);

			// Load initial electron temperature profile.
			real_t *T_init = LoadDataR(
				modulename,
				eqsys->GetFluidGrid()->GetRadialGrid(),
				s, "init"
			);

			if (!eqsys->IsAssigningToAlternativeEquation(id_T_cold)) {
				// NOTE: We should still mark the initial value as used,
				// to prevent complaints by the Python settings interface
				// when loading back the settings from the output file.
				if (T_init == nullptr)
					throw SettingsException("No initial data loaded for T_cold (from '%s/init'). Perhaps it has not been provided correctly?", modulename.c_str());
				eqsys->SetInitialValue(id_T_cold, T_init);
			}

			if (T_init != nullptr)
				delete [] T_init;
		} break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_T_COLD, type 
        );
    }
}

/**
 * Construct the equation for a prescribed temperature.
 */
void SimulationGenerator::ConstructEquation_T_cold_prescribed(
    const string& modulename, EquationSystem *eqsys, Settings *s,
	const len_t id_T
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Operator *eqn = new FVM::Operator(fluidGrid);

    FVM::Interpolator1D *interp = LoadDataRT_intp(modulename,fluidGrid->GetRadialGrid(), s);
    eqn->AddTerm(new FVM::PrescribedParameter(fluidGrid, interp));

    eqsys->SetOperator(id_T, id_T, eqn, "Prescribed");

    // Initialization
	eqsys->SetInitialValue(id_T, interp->Eval(0));
}


/**
 * Construct the equation for a self-consistent temperature evolution.
 */
void SimulationGenerator::ConstructEquation_T_cold_selfconsistent(
	const string& modulename,
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
    struct OtherQuantityHandler::T_terms *oqty_terms,
	const len_t id_eqn, const len_t id_T, const len_t id_W, const len_t id_n,
	const len_t id_j, bool isForThot
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    IonHandler *ionHandler = eqsys->GetIonHandler();

    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    /*len_t id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    len_t id_W_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_n_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);*/

    len_t id_E_field = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);


    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);

    Op1->AddTerm(new FVM::TransientTerm(fluidGrid,id_W));

    // Check if halo region heat losses should be included
    bool lcfs_user_input_psi = (len_t)s->GetInteger(MODULENAME_NRE  "/lcfs_user_input_psi");
    real_t lcfs_psi_edge_t0 = s->GetReal(MODULENAME_NRE "/lcfs_psi_edge_t0");

    bool parallel_losses = s->GetBool(modulename + "/halo_region_losses");
    if (parallel_losses) {
        HaloRegionHeatLossTerm* Par = new HaloRegionHeatLossTerm(fluidGrid,unknowns,ionHandler,-1, id_T, lcfs_user_input_psi, lcfs_psi_edge_t0);
        oqty_terms->halo = Par;
        Op1->AddTerm(Par); // Add the term for parallel losses
    }

    oqty_terms->ohmic = new OhmicHeatingTerm(fluidGrid, unknowns, id_j);
    Op2->AddTerm(oqty_terms->ohmic);

    bool withRecombinationRadiation = s->GetBool(modulename + "/recombination");

    // Load opacity settings
    len_t nitypes;
    const int_t *iopacity_modes = s->GetIntegerArray(MODULENAME_ION "/opacity_modes", 1, &nitypes);
    enum OptionConstants::ion_opacity_mode *opacity_mode = new enum OptionConstants::ion_opacity_mode[nitypes];
    for (len_t i = 0; i < nitypes; i++)
        opacity_mode[i] = (enum OptionConstants::ion_opacity_mode)iopacity_modes[i];

    oqty_terms->radiation = new RadiatedPowerTerm(
        fluidGrid,unknowns,ionHandler,adas,nist,amjuel,
		id_T, id_n,
		opacity_mode,withRecombinationRadiation
    );
    Op3->AddTerm(oqty_terms->radiation);


    FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
    // Add transport terms, if enabled
    bool hasTransport = ConstructTransportTerm(
        Op4, modulename, fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI,
        eqsys, s, false, true,
        &oqty_terms->advective_bc,
		&oqty_terms->diffusive_bc,
		nullptr, "transport", id_T, id_n
    );

    eqsys->SetOperator(id_eqn, id_E_field,Op2);
    eqsys->SetOperator(id_eqn, id_n,Op3);
    string desc = "dWc/dt = j_ohm*E - sum_i n_cold*n_i*L_i";

    bool includeNBI = false;
    if (s->HasSetting(MODULENAME "/NBI/enabled"))
        includeNBI = s->GetBool(MODULENAME "/NBI/enabled");

    if (includeNBI) {
        real_t s_max = s->GetReal(MODULENAME "/NBI/s_max");
        real_t r_beam = s->GetReal(MODULENAME "/NBI/r_beam");
        real_t Ti_beam = s->GetReal(MODULENAME "/NBI/Ti_beam");
        real_t m_i_beam = s->GetReal(MODULENAME "/NBI/m_i_beam");
        real_t R0 = s->GetReal(MODULENAME "/NBI/R0");
        len_t dims[1];
        const real_t *P0 = s->GetRealArray(MODULENAME "/NBI/P0", 1, dims);
        const real_t *n = s->GetRealArray(MODULENAME "/NBI/n", 1, dims);
        const real_t *energy_fractions = s->GetRealArray(MODULENAME "/NBI/energy_fractions", 1, dims);
        int_t gaussian_profile = s->GetInteger(MODULENAME "/NBI/gaussian_profile");

        // Load time-dependent j_B data
        FVM::Interpolator1D *j_B_profile = LoadDataT(
            modulename + "/NBI", s, "j_B"
        );

        FVM::Interpolator1D *Power_Profile = LoadDataT(
            modulename + "/NBI", s, "P_NBI"
        );

        NBIHandler *handler = eqsys->NBI_handler;
        if (handler == nullptr) {
            handler = new NBIHandler(eqsys->GetFluidGrid(), adas, ionHandler);
            handler->ConfigureFromSettings(
                s, unknowns,
                s_max, r_beam,
                P0, n, energy_fractions,
                Ti_beam, m_i_beam,
                j_B_profile, R0,
                gaussian_profile, Power_Profile
            );
            
        }
        eqsys->NBI_handler = handler;
  
        auto *nbi_e = new NBIElectronTerm(handler, eqsys->GetFluidGrid(), unknowns, ionHandler);
        Op4->AddTerm(nbi_e);
        oqty_terms->NBI = nbi_e;
        desc += " + NBI";
        eqsys->SetOperator(id_eqn, id_T, Op4); 
    }

    // SPI heat absorbtion
    OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode = (enum OptionConstants::eqterm_spi_heat_absorbtion_mode)s->GetInteger(MODULENAME_SPI "/heatAbsorbtion");
    if(spi_heat_absorbtion_mode!=OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT){
        Op4->AddTerm(new SPIHeatAbsorbtionTerm(fluidGrid,eqsys->GetSPIHandler(),-1));
    }

    // SPI ionization losses (from neutral to equilibrium)

    // TODO: simplify the bool logic below if possible
    OptionConstants::eqterm_ionization_mode ionization_mode =
        (enum OptionConstants::eqterm_ionization_mode)s->GetInteger(MODULENAME_ION "/ionization");
    bool includeKineticIonization = (ionization_mode == OptionConstants::EQTERM_IONIZATION_MODE_KINETIC)
                                    || (ionization_mode == OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC);
    if(includeKineticIonization && !(eqsys->HasHotTailGrid()||eqsys->HasRunawayGrid()))
        throw SettingsException("Invalid ionization mode: cannot use kinetic ionization without a kinetic grid.");
    bool collfreqModeIsFull = (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode")
    == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL;
    bool addFluidIonization = !(includeKineticIonization && eqsys->HasHotTailGrid() && collfreqModeIsFull);
    bool addFluidJacobian = (includeKineticIonization && eqsys->HasHotTailGrid() && (ionization_mode==OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC));

    len_t nZSPInShard;
    OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode = (enum OptionConstants::eqterm_spi_deposition_mode)s->GetInteger(MODULENAME_SPI "/deposition");
    OptionConstants::eqterm_spi_abl_ioniz_mode spi_abl_ioniz_mode = (enum OptionConstants::eqterm_spi_abl_ioniz_mode)s->GetInteger(MODULENAME_SPI "/abl_ioniz");
    const real_t *SPIMolarFraction  = s->GetRealArray(MODULENAME_ION "/SPIMolarFraction", 1, &nZSPInShard);

    // Add one ionization loss term for every species the pellet consists of
    // Note that, when accounting for the heat absorbed in the neutral cloud,
    // this energy is currently being redeposited where the material is deposited,
    // and we therefore need these terms even when heat absorption is included.
    if(spi_deposition_mode!=OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT && spi_abl_ioniz_mode!=OptionConstants::EQTERM_SPI_ABL_IONIZ_MODE_NEUTRAL){
        len_t offset=0;
        len_t nShard = eqsys->GetSPIHandler()->GetNShard();
        const len_t nZ = ionHandler->GetNZ();
        for(len_t iZ=0;iZ<nZ;iZ++){
            if(SPIMolarFraction[offset]>0){
                Op4->AddTerm(new IonSPIIonizLossTerm(fluidGrid, eqsys->GetIonHandler(), iZ, adas, eqsys->GetUnknownHandler(),
                                                     addFluidIonization, addFluidJacobian, eqsys->GetSPIHandler(), SPIMolarFraction,offset,1,nist,false, spi_abl_ioniz_mode));
                offset+=nShard;
            }else {
                offset+=1;
            }
        }
    }

    if(hasTransport){
        oqty_terms->transport = Op4->GetAdvectionDiffusion();
        desc += " + transport";
    }

    if (spi_deposition_mode!=OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT ||
        spi_heat_absorbtion_mode!=OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT ||
        hasTransport)
        eqsys->SetOperator(id_eqn, id_T,Op4);

    /**
     * If hot-tail grid is enabled, add collisional
     * energy transfer from hot-tail to T_cold.
     * NOTE: We temporarily disable the collisional energy transfer
     *       in collfreq_mode FULL because we lack the correction
     *       term for when electrons transfer from the cold to the
     *       hot region. This should be corrected for at some point!
     */
	// When building the equation for T_hot, we want to disable this term
	// altogether since the Maxwellian should not heat itself.
    bool collfreqModeFull = ((enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode") == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);
    if (eqsys->HasHotTailGrid()&& !collfreqModeFull && !isForThot) {
        len_t id_f_hot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);

        FVM::MomentQuantity::pThresholdMode pMode =
            (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
        real_t pThreshold = 0.0;
        if(collfreqModeFull){
            // With collfreq_mode FULL, only add contribution from hot electrons
            // defined as those with momentum above the defined threshold.
            pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        }
        oqty_terms->fhot_coll = new CollisionalEnergyTransferKineticTerm(
            fluidGrid,eqsys->GetHotTailGrid(),
            id_eqn, id_f_hot,eqsys->GetHotTailCollisionHandler(), eqsys->GetUnknownHandler(),
            eqsys->GetHotTailGridType(), -1.0,
            pThreshold, pMode
        );
        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        Op5->AddTerm(oqty_terms->fhot_coll );
        eqsys->SetOperator(id_eqn, id_f_hot, Op5);
        desc += " + int(W*nu_E*f_hot)";
    }

	// Neglect RE heating of T_hot
	if (!isForThot) {
		// If runaway grid and not FULL collfreqmode, add collisional
		// energy transfer from runaways to T_cold.
		if (eqsys->HasRunawayGrid()) {
			len_t id_f_re = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);

			oqty_terms->fre_coll = new CollisionalEnergyTransferKineticTerm(
				fluidGrid,eqsys->GetRunawayGrid(),
				id_eqn, id_f_re,eqsys->GetRunawayCollisionHandler(),eqsys->GetUnknownHandler(),
				eqsys->GetRunawayGridType(), -1.0
			);
			FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
			Op5->AddTerm(oqty_terms->fre_coll );
			eqsys->SetOperator(id_eqn, id_f_re, Op5);
			desc += " + int(W*nu_E*f_re)";
		} else {
			len_t id_n_re = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
			oqty_terms->nre_coll = new CollisionalEnergyTransferREFluidTerm(
				fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid()->GetLnLambda(), -1.0
			);
			FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
			Op5->AddTerm(oqty_terms->nre_coll );
			eqsys->SetOperator(id_eqn, id_n_re, Op5);

			desc += " + e*c*Ec*n_re";
		}
	}

    // ADD COLLISIONAL ENERGY TRANSFER WITH ION SPECIES
    OptionConstants::uqty_T_i_eqn Ti_type =
        (OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi");
    if(Ti_type == OptionConstants::UQTY_T_I_INCLUDE) {
        CoulombLogarithm *lnLambda;
		if (isForThot)
			lnLambda = eqsys->GetREFluid()->GetLnLambdaHot();
		else
			lnLambda = eqsys->GetREFluid()->GetLnLambda();

        const len_t nZ = ionHandler->GetNZ();
        const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER);
        oqty_terms->ion_coll = new FVM::Operator(fluidGrid);
        for(len_t iz=0; iz<nZ; iz++){
            oqty_terms->ion_coll->AddTerm(
                new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    0, false,
                    iz, true,
					id_T, id_W, id_n,
                    unknowns, lnLambda, ionHandler, -1.0)
                );
        }
        eqsys->SetOperator(id_eqn, id_Wi, oqty_terms->ion_coll);
        desc += " + sum_i Q_ei";
    }
    eqsys->SetOperator(id_eqn, id_W,Op1,desc);
}


/**
 * Construct the equation for electron energy content:
 *    W_cold = 3n_cold*T_cold/2
 */
void SimulationGenerator::ConstructEquation_W_cold(
    EquationSystem *eqsys, Settings* /*s*/
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);

    len_t id_W_cold = eqsys->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_T_cold = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);
    len_t id_n_cold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);

    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0) );
    Op2->AddTerm(new ElectronHeatTerm(fluidGrid,id_n_cold,eqsys->GetUnknownHandler()) );

    eqsys->SetOperator(id_W_cold, id_W_cold, Op1, "W_cold = (3/2)*n_cold*T_cold");
    eqsys->SetOperator(id_W_cold, id_T_cold, Op2);

    eqsys->initializer->AddRule(
        id_W_cold,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_T_cold,
        id_n_cold 
    );
}
