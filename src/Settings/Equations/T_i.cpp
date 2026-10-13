#include "DREAM/Equations/Fluid/IonSpeciesIdentityTerm.hpp"
#include "DREAM/Equations/Fluid/IonSpeciesTransientTerm.hpp"
#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Operator.hpp"
#include "DREAM/Equations/Fluid/NBIIonTerm.hpp"
#include "DREAM/NBIHandler.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/Fluid/NeutralCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Equations/Fluid/ExotermicEnergyTransferTerm.hpp"

/**
 * Implementation of equations governing the evolution of the
 * ion heat W_i = 1.5 * N_i * T_i, which satisfy
 *   dW_i/dt = Q_ie + sum_j Q_ij
 * with the collisional energy transfer Q_ij summed over all ion species j
 * and the exchange Q_ie with electrons. This ion heat is not
 * charge-state resolved, so that each species only has one temperature.
 */

using namespace DREAM;
using namespace std;


#define MODULENAME "eqsys/n_i"
#define MODULENAME_TCOLD "eqsys/T_cold"


void SimulationGenerator::ConstructEquation_T_i(EquationSystem *eqsys, Settings *s, struct OtherQuantityHandler::eqn_terms *oqty_terms){
    /**
     * if the electron heat W_cold is evolved self-consistently,
     * also evolve the ion heat W_i. Otherwise set it to constant.
     */
    enum OptionConstants::uqty_T_cold_eqn TcoldType = (enum OptionConstants::uqty_T_cold_eqn)s->GetInteger("eqsys/T_cold/type");
    if(TcoldType==OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED)
        ConstructEquation_T_i_trivial(eqsys, s);
    else if (TcoldType == OptionConstants::UQTY_T_COLD_SELF_CONSISTENT)
        ConstructEquation_T_i_selfconsistent(eqsys, s, oqty_terms);
    else 
        throw SettingsException(
            "T_i: Unrecognized equation type for '%s': %d.",
            OptionConstants::UQTY_T_COLD, TcoldType
        );

    // Initialize heat from ion densities and input ion temperatures
    real_t *Ti_data = LoadDataIonR(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s, eqsys->GetIonHandler()->GetNZ(), "initialTi");

    const len_t nValues = eqsys->GetIonHandler()->GetNZ()* eqsys->GetFluidGrid()->GetNr();
    const std::vector<real_t> Ti_init(Ti_data, Ti_data + nValues);
    delete [] Ti_data;

	len_t id_Ni = eqsys->GetUnknownID(OptionConstants::UQTY_NI_DENS);
	len_t id_ions = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
	const bool separateNeutrals = s->GetInteger("eqsys/n_i/typeTn") == OptionConstants::UQTY_T_N_INCLUDE;

	std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_Ti =
		[Ti_init,id_Ni,id_ions,eqsys,separateNeutrals](FVM::UnknownQuantityHandler *u, real_t *Ti) {
		const real_t *Ni_init = u->GetUnknownInitialData(id_Ni);
		const real_t *ions_init = u->GetUnknownInitialData(id_ions);
		const len_t nZ = eqsys->GetIonHandler()->GetNZ();
		const len_t nr = eqsys->GetFluidGrid()->GetNr();
		for(len_t iz=0; iz<nZ; iz++) {
			for(len_t ir=0; ir<nr; ir++) {
				const len_t it = iz*nr + ir;
				real_t n = Ni_init[it];
                //If we want to separate the neutrals, 
                //we need to sum over all other charge states to get the 
                //total density of the species (remove neutral density).
				if (separateNeutrals) {
					n = 0;
					for (len_t Z0=1; Z0<=eqsys->GetIonHandler()->GetZ(iz); Z0++) {
						const len_t indZ = eqsys->GetIonHandler()->GetIndex(iz, Z0);
						n += ions_init[indZ*nr + ir];
					}
				}
				Ti[it] = Ti_init[it] * 1.5*Constants::ec*n;
			}
		}
	};
    //eqsys->SetInitialValue(eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER), Ti_init);    
	len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER);
	eqsys->initializer->AddRule(
		id_Wi,
		EqsysInitializer::INITRULE_EVAL_FUNCTION,
		initfunc_Ti,
		id_Ni, id_ions
	);
}

/**
 * Set a trivial equation W_i = constant for all ion species
 */
void SimulationGenerator::ConstructEquation_T_i_trivial(EquationSystem *eqsys, Settings* /*s*/){
    const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER); 
    FVM::Operator *Op = new FVM::Operator(eqsys->GetFluidGrid());
    for(len_t iz=0; iz<eqsys->GetIonHandler()->GetNZ(); iz++)
        Op->AddTerm( 
            new IonSpeciesTransientTerm(eqsys->GetFluidGrid(), iz, id_Wi)
        );
    eqsys->SetOperator(id_Wi, id_Wi, Op, "dW_i/dt = 0");
}


/** 
 * Implements the self-consistent evolution of ion heat W_i for each species
 */
void SimulationGenerator::ConstructEquation_T_i_selfconsistent(EquationSystem *eqsys, Settings* s,  struct OtherQuantityHandler::eqn_terms *oqty_terms){
    const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER); 
    const len_t id_Wcold = eqsys->GetUnknownID(OptionConstants::UQTY_W_COLD);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    IonHandler *ionHandler = eqsys->GetIonHandler();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();    
    const len_t nZ = ionHandler->GetNZ();

    const bool separateNeutrals = s->GetInteger("eqsys/n_i/typeTn") == OptionConstants::UQTY_T_N_INCLUDE;

    FVM::Operator *Op_Wij = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_Wie = new FVM::Operator(fluidGrid);

    CoulombLogarithm *lnLambda = eqsys->GetREFluid()->GetLnLambda();
    NBIHandler *handler = eqsys->NBI_handler;

    
    oqty_terms->T_i_Qij.resize(nZ);
    for(len_t iz=0; iz<nZ; iz++)
        oqty_terms->T_i_Qij[iz].resize(nZ, nullptr);
    oqty_terms->T_i_Qie.resize(nZ, nullptr);


    for(len_t iz=0; iz<nZ; iz++){
        Op_Wij->AddTerm(
            new IonSpeciesTransientTerm(fluidGrid, iz, id_Wi, -1.0)
        );
        for(len_t jz=0; jz<nZ; jz++){
            if(jz==iz) // the term is trivial =0 for self collisions and can be skipped
                continue;
            auto *Qij = new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    iz, true,
                    jz, true,
                    unknowns, lnLambda, ionHandler, 1.0, separateNeutrals
                );
            Op_Wij->AddTerm(Qij);
            oqty_terms->T_i_Qij[iz][jz] = Qij;
        }
        auto *Qie = new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    iz, true,
                    0, false,
                    unknowns, lnLambda, ionHandler, 1.0, separateNeutrals
                );
        Op_Wie->AddTerm(Qie);
        oqty_terms->T_i_Qie[iz] = Qie;
    }

    bool includeNBI = false;
    if (s->HasSetting(MODULENAME_TCOLD "/NBI/enabled")) {
        includeNBI = s->GetBool(MODULENAME_TCOLD "/NBI/enabled");
    }
    if (includeNBI){
        oqty_terms->T_i_NBI.resize(nZ, nullptr);
        for(len_t iz=0; iz<nZ; iz++){
           auto *nbi_i = new NBIIonTerm(handler, fluidGrid, ionHandler, unknowns, iz);
           Op_Wij->AddTerm(nbi_i);
           oqty_terms->T_i_NBI[iz] = nbi_i;           
       }
    }

    eqsys->SetOperator(id_Wi, id_Wi, Op_Wij, "dW_i/dt = sum_j Q_ij + Q_ie");
    eqsys->SetOperator(id_Wi, id_Wcold, Op_Wie);
}

/*
 * Called from EquationSystem.cpp to construct the equation for the neutral heat W_n if enabled.
 */
void SimulationGenerator::ConstructEquation_T_n(
    EquationSystem *eqsys,
    Settings *s,
    struct OtherQuantityHandler::eqn_terms*
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    IonHandler *ions = eqsys->GetIonHandler();

    const len_t id_Wn = eqsys->GetUnknownID( OptionConstants::UQTY_WN_ENER);

    const len_t nSpecies = ions->GetNZ();

    FVM::Operator *Op_Wn = new FVM::Operator(fluidGrid);

    for (len_t iz = 0; iz < nSpecies; iz++) {
       Op_Wn->AddTerm(
          new IonSpeciesTransientTerm(
              fluidGrid, iz, id_Wn, -1.0
          )
      );
        for (len_t jz = iz + 1; jz < nSpecies; jz++) {
            Op_Wn->AddTerm(
                new NeutralCollisionalEnergyTransferTerm(
                    fluidGrid, iz, jz, unknowns, ions
                )
            );
        }
        Op_Wn->AddTerm(
              new ExotermicEnergyTransferTerm(
                  fluidGrid, iz, unknowns, ions, eqsys->GetRateHandler()
              )
          );
  }

    eqsys->SetOperator(id_Wn, id_Wn, Op_Wn, "dW_n/dt = neutral thermal terms");

     // Load the initial neutral temperature, in eV.
    real_t *Tn_data = LoadDataIonR("eqsys/n_i", fluidGrid->GetRadialGrid(),s, nSpecies,"initialTn");

    if (Tn_data == nullptr) {
        throw SettingsException(
            "No initial neutral temperature loaded."
        );
    }

    const len_t nr = fluidGrid->GetNr();
    const len_t id_ions = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);


  // Keep the temperatures alive until the initializer runs.
    const std::vector<real_t> Tn_init(Tn_data, Tn_data + nSpecies * nr);
    delete [] Tn_data;


    // Calculates only neutral thermal energy
    std::function<void(FVM::UnknownQuantityHandler*, real_t*)>
    initializeNeutralEnergy =
    [Tn_init, id_ions, ions, nSpecies, nr](
        FVM::UnknownQuantityHandler *unknowns,
          real_t *Wn
    ) {
        const real_t *initialDensities =
            unknowns->GetUnknownInitialData(id_ions);

        for (len_t iz = 0; iz < nSpecies; iz++) {
            const len_t neutralIndex = ions->GetIndex(iz, 0);

            for (len_t ir = 0; ir < nr; ir++) {
                const len_t energyIndex = iz * nr + ir;
                const len_t densityIndex = neutralIndex * nr + ir;

                const real_t neutralDensity = initialDensities[densityIndex];

                Wn[energyIndex] =1.5 * Constants::ec * neutralDensity * Tn_init[energyIndex];
              }
          }
      };

  eqsys->initializer->AddRule(
      id_Wn,
      EqsysInitializer::INITRULE_EVAL_FUNCTION,
      initializeNeutralEnergy,
      id_ions
  );
  }
