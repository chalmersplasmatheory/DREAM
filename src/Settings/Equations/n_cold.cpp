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

#define MODULENAME "eqsys/n_cold"

/**
 * Construct the equation for the cold electron density, 'n_cold'.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct
 *        the equation.
 */
void SimulationGenerator::ConstructEquation_n_cold(
    EquationSystem *eqsys, Settings *s
) {
    enum OptionConstants::uqty_n_cold_eqn eqn = (enum OptionConstants::uqty_n_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (eqn) {
        case OptionConstants::UQTY_N_COLD_EQN_PRESCRIBED:
            ConstructEquation_n_cold_prescribed(eqsys, s);
            break;
        
        case OptionConstants::UQTY_N_COLD_EQN_SELFCONSISTENT:
            ConstructEquation_n_cold_selfconsistent(eqsys, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for 'n_cold': %d.",
                eqn
            );
    }
}

/**
 * Construct the equation describing a prescribed cold
 * electron density.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct
 *        the equation.
 */
void SimulationGenerator::ConstructEquation_n_cold_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    // const real_t t0 = 0;
    FVM::Operator *eqn = new FVM::Operator(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT_intp(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetOperator(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_COLD, eqn, "Prescribed");
    //eqsys->SetInitialValue(OptionConstants::UQTY_N_COLD, interp->Eval(t0), t0);

    eqsys->initializer->AddRule(
        OptionConstants::UQTY_N_COLD,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );
}

/**
 * Construct the equation describing the cold electron density as
 * the free electron density that is not hot or runaway
 */
void SimulationGenerator::ConstructEquation_n_cold_selfconsistent(
    EquationSystem *eqsys, Settings */*s*/
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    const len_t id_ncold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const len_t id_nhot  = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
    const len_t id_nre   = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    const len_t id_ni    = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);


    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
    
    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid));
    Op2->AddTerm(new FVM::IdentityTerm(fluidGrid));
    Op3->AddTerm(new FVM::IdentityTerm(fluidGrid));
    Op4->AddTerm(new FreeElectronDensityTerm(fluidGrid,eqsys->GetIonHandler(),-1.0));
    eqsys->SetOperator(id_ncold, id_ncold, Op1, "n_cold = n_free - n_hot - n_re");
    eqsys->SetOperator(id_ncold, id_nhot, Op2);
    eqsys->SetOperator(id_ncold, id_nre, Op3);
    eqsys->SetOperator(id_ncold, id_ni, Op4);

    // Initialization
    eqsys->initializer->AddRule(
        OptionConstants::UQTY_N_COLD,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_ni, id_nre, id_nhot
    );



}

