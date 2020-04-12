/**
 * Definition of equations relating to n_cold.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;

#define MODULENAME "equationsystem/n_cold"

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
    const real_t t0 = 0;
    FVM::Equation *eqn = new FVM::Equation(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetEquation(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_COLD, eqn);
    //eqsys->SetInitialValue(OptionConstants::UQTY_N_COLD, interp->Eval(t0), t0);
}

