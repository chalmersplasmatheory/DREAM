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
    const real_t *Ni_init = eqsys->GetUnknownHandler()->GetUnknownInitialData(eqsys->GetUnknownID(OptionConstants::UQTY_NI_DENS));
    for(len_t it=0; it<eqsys->GetIonHandler()->GetNZ()*eqsys->GetFluidGrid()->GetNr(); it++)
        Ti_init[it] *= 1.5*Constants::ec*Ni_init[it];
    eqsys->SetInitialValue(eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER), Ti_init);    
}

/**
 * Set a trivial equation W_i = constant for all ion species
 */
void SimulationGenerator::ConstructEquation_T_i_trivial(EquationSystem *eqsys, Settings */*s*/){
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
void SimulationGenerator::ConstructEquation_T_i_selfconsistent(EquationSystem *eqsys, Settings */*s*/){
    const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER); 
    const len_t id_Ni = eqsys->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    const len_t id_ncold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
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
        const len_t Zi  = ionHandler->GetZ(iz);
        const real_t mi = ionHandler->GetIonSpeciesMass(iz);
        for(len_t jz=0; jz<nZ; jz++){
            if(jz==iz) // the term is trivial =0 for self collisions and can be skipped
                continue;
            const len_t Zj  = ionHandler->GetZ(jz);
            const real_t mj = ionHandler->GetIonSpeciesMass(jz);
            Op_Wij->AddTerm(
                new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid, 
                    id_Ni, id_Wi, Zi, mi, iz, 
                    id_Ni, id_Wi, Zj, mj, jz,
                    unknowns, lnLambda, false)
            );
        }
        Op_Wie->AddTerm(
            new MaxwellianCollisionalEnergyTransferTerm(
                fluidGrid, 
                id_Ni, id_Wi, Zi, mi, iz, 
                id_ncold, id_Wcold, 1, Constants::me, 0,
                unknowns, lnLambda, true)
        );
    }
    eqsys->SetOperator(id_Wi, id_Wi, Op_Wij, "dW_i/dt = sum_j Q_ij + Q_ie");
    eqsys->SetOperator(id_Wi, id_Wcold, Op_Wie);
}
