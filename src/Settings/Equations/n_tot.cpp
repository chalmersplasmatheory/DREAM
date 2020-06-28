/**
 * Set n_tot (total electron density) from quasi-neutrality.
 */

#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/TotalElectronDensityTerm.hpp"
#include "FVM/Equation/IdentityTerm.hpp"


using namespace DREAM;


/**
 * Construct the equation evaluating n_tot.
 */
void SimulationGenerator::ConstructEquation_n_tot(
    EquationSystem *eqsys, Settings*
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);

    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
    Op2->AddTerm(new TotalElectronDensityTerm(fluidGrid,eqsys->GetIonHandler()));
    
    eqsys->SetOperator(OptionConstants::UQTY_N_TOT, OptionConstants::UQTY_N_TOT, Op1, "n_tot = sum_i Z_i n_i");
    eqsys->SetOperator(OptionConstants::UQTY_N_TOT, OptionConstants::UQTY_ION_SPECIES, Op2);
    
    // Initialization
    const len_t id_n_i = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    eqsys->initializer->AddRule(
        OptionConstants::UQTY_N_TOT,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_n_i
    );
}

