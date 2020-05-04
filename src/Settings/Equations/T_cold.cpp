/**
 * Definition of equations relating to the cold electron temperature.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


#define MODULENAME "equationsystem/T_cold"

/**
 * Construct the equation for the electric field.
 */
void SimulationGenerator::ConstructEquation_T_cold(
    EquationSystem *eqsys, Settings *s
) {
    enum OptionConstants::uqty_T_cold_eqn type = (enum OptionConstants::uqty_T_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED:
            ConstructEquation_T_cold_prescribed(eqsys, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_T_COLD, type
            );
    }
}

/**
 * Construct the equation for a prescribed electric field.
 */
void SimulationGenerator::ConstructEquation_T_cold_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    // const real_t t0 = 0;
    FVM::Equation *eqn = new FVM::Equation(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetEquation(OptionConstants::UQTY_T_COLD, OptionConstants::UQTY_T_COLD, eqn);
}

