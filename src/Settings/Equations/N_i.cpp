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
 * simply summing the densities of all charge states
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
    
    // Operators for net ion density equation
    FVM::Operator *Op_Ni = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_ni = new FVM::Operator(fluidGrid);    
    for(len_t i=0; i<nZ; i++)
        switch(types[i]) {
            case OptionConstants::ION_DATA_PRESCRIBED:
                [[fallthrough]]
            case OptionConstants::ION_DATA_TYPE_DYNAMIC:
                Op_Ni->AddTerm( 
                    new IonSpeciesIdentityTerm(fluidGrid, i, -1.0) 
                );
                Op_ni->AddTerm( 
                    new NetIonDensityFromIonChargeStatesTerm(fluidGrid, ih->GetZ(i), i, ih) 
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
    len_t nr = fluidGrid->GetNr();
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

    eqsys->SetInitialValue(id_Ni, Ni0);
    delete [] types;
}