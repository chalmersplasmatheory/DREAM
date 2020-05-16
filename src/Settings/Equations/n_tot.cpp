/**
 * Set n_tot (total electron density) from quasi-neutrality.
 */

#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/NTotFromQuasiNeutrality.hpp"
#include "FVM/Equation/IdentityTerm.hpp"


using namespace DREAM;


/**
 * Construct the equation evaluating n_tot.
 */
void SimulationGenerator::ConstructEquation_n_tot(
    EquationSystem *eqsys, Settings*
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Equation *eqn = new FVM::Equation(fluidGrid);

    eqn->AddTerm(new NTotFromQuasiNeutrality(fluidGrid, eqsys->GetIonHandler()));
    eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

    eqsys->SetEquation(OptionConstants::UQTY_N_TOT, OptionConstants::UQTY_N_TOT, eqn);
}

