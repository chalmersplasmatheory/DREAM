/**
 * Definition of equations relating to n_cold.
 */

#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/FreeElectronDensityTerm.hpp"
#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;


/**
 * Construct the equation describing the cold electron density as
 * the free electron density that is not hot or runaway
 */
void SimulationGenerator::ConstructEquation_n_abl(
    EquationSystem *eqsys, Settings */*s*/
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    const len_t id_nabl = eqsys->GetUnknownID(OptionConstants::UQTY_N_ABL);
    const len_t id_ni_abl    = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES_ABL);

	FVM::Operator *Opid = new FVM::Operator(fluidGrid);
    FVM::Operator *Op = new FVM::Operator(fluidGrid);
    
    Opid->AddTerm(new FVM::IdentityTerm(fluidGrid));
    Op->AddTerm(new FreeElectronDensityTerm(fluidGrid,eqsys->GetIonHandler(),-1.0));
    eqsys->SetOperator(id_nabl, id_nabl, Opid, "dn_abl/dt = ionisation from SPI material ");
    eqsys->SetOperator(id_nabl, id_ni_abl, Op);

    // Initialization
    eqsys->initializer->AddRule(
        OptionConstants::UQTY_N_ABL,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_ni_abl
    );



}

