#include "DREAM/Equations/Fluid/IonSpeciesIdentityTerm.hpp"
#include "DREAM/Equations/Fluid/IonSpeciesTransientTerm.hpp"
#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Operator.hpp"

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


void SimulationGenerator::ConstructEquation_T_i(EquationSystem *eqsys, Settings *s){
    /**
     * if the electron heat W_cold is evolved self-consistently,
     * also evolve the ion heat W_i. Otherwise set it to constant.
     */
    enum OptionConstants::uqty_T_cold_eqn TcoldType = (enum OptionConstants::uqty_T_cold_eqn)s->GetInteger("eqsys/T_cold/type");
    if(TcoldType==OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED)
        ConstructEquation_T_i_trivial(eqsys, s);
    else if (TcoldType == OptionConstants::UQTY_T_COLD_SELF_CONSISTENT)
        ConstructEquation_T_i_selfconsistent(eqsys, s);
    else 
        throw SettingsException(
            "T_i: Unrecognized equation type for '%s': %d.",
            OptionConstants::UQTY_T_COLD, TcoldType
        );

    // Initialize heat from ion densities and input ion temperatures
    real_t *Ti_init = LoadDataIonR(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s, eqsys->GetIonHandler()->GetNZ(), "initialTi");
	len_t id_Ni = eqsys->GetUnknownID(OptionConstants::UQTY_NI_DENS);

	std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_Ti =
		[Ti_init,id_Ni,eqsys](FVM::UnknownQuantityHandler *u, real_t *Ti) {
		const real_t *Ni_init = u->GetUnknownInitialData(id_Ni);
		for(len_t it=0; it<eqsys->GetIonHandler()->GetNZ()*eqsys->GetFluidGrid()->GetNr(); it++)
			Ti[it] = Ti_init[it] * 1.5*Constants::ec*Ni_init[it];
	};
    //eqsys->SetInitialValue(eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER), Ti_init);    
	len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER);
	eqsys->initializer->AddRule(
		id_Wi,
		EqsysInitializer::INITRULE_EVAL_FUNCTION,
		initfunc_Ti,
		id_Ni
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
void SimulationGenerator::ConstructEquation_T_i_selfconsistent(EquationSystem *eqsys, Settings* /*s*/){
    const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER); 
    const len_t id_Wcold = eqsys->GetUnknownID(OptionConstants::UQTY_W_COLD);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    IonHandler *ionHandler = eqsys->GetIonHandler();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();    
    const len_t nZ = ionHandler->GetNZ();


    FVM::Operator *Op_Wij = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_Wie = new FVM::Operator(fluidGrid);

    CoulombLogarithm *lnLambda = eqsys->GetREFluid()->GetLnLambda();
    for(len_t iz=0; iz<nZ; iz++){
        Op_Wij->AddTerm( 
            new IonSpeciesTransientTerm(fluidGrid, iz, id_Wi, -1.0)
        );
        for(len_t jz=0; jz<nZ; jz++){
            if(jz==iz) // the term is trivial =0 for self collisions and can be skipped
                continue;
            Op_Wij->AddTerm(
                new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    iz, true,
                    jz, true,
                    unknowns, lnLambda, ionHandler)
            );
        }
        Op_Wie->AddTerm(
            new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    iz, true,
                    0, false,
                    unknowns, lnLambda, ionHandler)
        );
    }
    eqsys->SetOperator(id_Wi, id_Wi, Op_Wij, "dW_i/dt = sum_j Q_ij + Q_ie");
    eqsys->SetOperator(id_Wi, id_Wcold, Op_Wie);
}
