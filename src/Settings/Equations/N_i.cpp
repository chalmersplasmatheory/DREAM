#include "DREAM/Equations/Fluid/IonSpeciesIdentityTerm.hpp"
#include "DREAM/Equations/Fluid/IonSpeciesTransientTerm.hpp"
#include "DREAM/Equations/Fluid/NetIonDensityFromIonChargeStatesTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/Operator.hpp"

/**
 * Implementation of equations governing the evolution of the
 * net density of ion species, satisfying
 *      N_i = sum_j n_i^(j)
 * simply summing the densities of all charge states.
 * With a separate neutral heat equation (W_n), the neutrals are left out
 * of N_i (so that W_i = 1.5*N_i*T_i), and are instead given by
 *      N_n = n_i^(0)
 */

using namespace DREAM;
using namespace std;

/**
 * Construct the equation governing the evolution of the
 * ion densities for each charge state.
 */
void SimulationGenerator::ConstructEquation_Ion_Ni(EquationSystem *eqsys, Settings *s) {
    len_t nZ;
    const int_t *itypes = s->GetIntegerArray("eqsys/n_i/types", 1, &nZ);
    enum OptionConstants::ion_data_type *types = new enum OptionConstants::ion_data_type[nZ];
    for (len_t i = 0; i < nZ; i++)
        types[i] = (enum OptionConstants::ion_data_type)itypes[i];

    IonHandler *ih = eqsys->GetIonHandler();
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid(); 

    const len_t id_ni = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    const len_t id_Ni = eqsys->GetUnknownID(OptionConstants::UQTY_NI_DENS);

    // With separate neutrals, N_i only contains the charged states (Z0 >= 1)
    const bool separateNeutrals = s->GetInteger("eqsys/n_i/typeTn") == OptionConstants::UQTY_T_N_INCLUDE;
    const len_t Z0min = separateNeutrals ? 1 : 0;
    
    // Operators for net ion density equation
    FVM::Operator *Op_Ni = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_ni = new FVM::Operator(fluidGrid);    
    for(len_t i=0; i<nZ; i++)
        switch(types[i]) {
            case OptionConstants::ION_DATA_PRESCRIBED:
                [[fallthrough]];
            case OptionConstants::ION_DATA_TYPE_DYNAMIC:
                Op_Ni->AddTerm( 
                    new IonSpeciesIdentityTerm(fluidGrid, i, -1.0) 
                );
                Op_ni->AddTerm( 
                    new NetIonDensityFromIonChargeStatesTerm(fluidGrid, ih->GetZ(i), i, ih, 1.0, Z0min)
                );
            break;
            default:
                throw SettingsException(
                    "N_i: Currently only supports PRESCRIBED and DYNAMIC. Provided: %d.",
                    types[i]
                );
                break;
        }
    eqsys->SetOperator(id_Ni, id_Ni, Op_Ni, "N_i = sum_j n_i^(j)");
    eqsys->SetOperator(id_Ni, id_ni, Op_ni);   
    
    // PLACEHOLDER: Does not account for equilibrium: for equilibrium species 
    // we will load initial values for Ni from settings, and then n_i will be 
    // initialized from those, somehow (probably by evaluating the corresponding
    // equation for the steady state charge distribution).
    /*len_t nr = fluidGrid->GetNr();
    const real_t *ni = eqsys->GetUnknownHandler()->GetUnknownInitialData(id_ni);
    real_t *Ni0 = new real_t[nZ*nr];
    for(len_t iz=0; iz<nZ; iz++)
        for(len_t ir=0; ir<nr; ir++){
            Ni0[iz*nr+ir] = 0;
            for(len_t Z0=0; Z0<=ih->GetZ(iz); Z0++){
                len_t indZ = ih->GetIndex(iz,Z0);
                Ni0[iz*nr+ir] += ni[indZ*nr+ir];
            }
        }
	*/

	std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_Ni =
		[fluidGrid,id_ni,nZ,ih,Z0min](FVM::UnknownQuantityHandler *u, real_t *Ni) {
		len_t nr = fluidGrid->GetNr();
		const real_t *ni = u->GetUnknownInitialData(id_ni);
		for(len_t iz=0; iz<nZ; iz++)
			for(len_t ir=0; ir<nr; ir++){
				Ni[iz*nr+ir] = 0;
				for(len_t Z0=Z0min; Z0<=ih->GetZ(iz); Z0++){
					len_t indZ = ih->GetIndex(iz,Z0);
					Ni[iz*nr+ir] += ni[indZ*nr+ir];
				}
			}
	};

    //eqsys->SetInitialValue(id_Ni, Ni0);
	eqsys->initializer->AddRule(
		id_Ni,
		EqsysInitializer::INITRULE_EVAL_FUNCTION,
		initfunc_Ni,
		id_ni
	);
    delete [] types;
}


/**
 * Construct the equation for the neutral density of each species,
 *      N_n = n_i^(0)
 * used with a separate neutral heat equation, so that W_n = 1.5*N_n*T_n.
 * Same as N_i above, but summing only the neutral charge state.
 */
void SimulationGenerator::ConstructEquation_Ion_Nn(EquationSystem *eqsys, Settings*) {
    IonHandler *ih = eqsys->GetIonHandler();
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    const len_t nZ = ih->GetNZ();

    const len_t id_ni = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    const len_t id_Nn = eqsys->GetUnknownID(OptionConstants::UQTY_NN_DENS);

    // Operators for the neutral density equation
    FVM::Operator *Op_Nn = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_ni = new FVM::Operator(fluidGrid);
    for (len_t i = 0; i < nZ; i++) {
        Op_Nn->AddTerm(
            new IonSpeciesIdentityTerm(fluidGrid, i, -1.0)
        );
        // Sum from Z0 = 0 to Z0 = 0, i.e. only the neutral state
        Op_ni->AddTerm(
            new NetIonDensityFromIonChargeStatesTerm(fluidGrid, 0, i, ih, 1.0, 0)
        );
    }
    eqsys->SetOperator(id_Nn, id_Nn, Op_Nn, "N_n = n_i^(0)");
    eqsys->SetOperator(id_Nn, id_ni, Op_ni);

    // Initialize from the initial neutral densities
    std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_Nn =
        [fluidGrid,id_ni,nZ,ih](FVM::UnknownQuantityHandler *u, real_t *Nn) {
        const len_t nr = fluidGrid->GetNr();
        const real_t *ni = u->GetUnknownInitialData(id_ni);
        for (len_t iz = 0; iz < nZ; iz++)
            for (len_t ir = 0; ir < nr; ir++)
                Nn[iz*nr+ir] = ni[ih->GetIndex(iz,0)*nr+ir];
    };

    eqsys->initializer->AddRule(
        id_Nn,
        EqsysInitializer::INITRULE_EVAL_FUNCTION,
        initfunc_Nn,
        id_ni
    );
}
