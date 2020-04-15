/**
 * Definition of equations relating to the electric field.
 *
 * Note that we define the electric field as
 *
 *      <E*B>
 *   -----------
 *   sqrt(<B^2>)
 *
 * in DREAM, where 'E' denotes the local electric field, 'B' the
 * local magnetic field and '<X>' denotes the flux-surface average
 * of a quantity 'X'.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


#define MODULENAME "equationsystem/E_field"

/**
 * Construct the equation for the electric field.
 */
void SimulationGenerator::ConstructEquation_E_field(
    EquationSystem *eqsys, Settings *s
) {
    enum OptionConstants::uqty_E_field_eqn type = (enum OptionConstants::uqty_E_field_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED:
            ConstructEquation_E_field_prescribed(eqsys, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_E_FIELD, type
            );
    }
}

/**
 * Construct the equation for a prescribed electric field.
 */
void SimulationGenerator::ConstructEquation_E_field_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    // const real_t t0 = 0;
    FVM::Equation *eqn = new FVM::Equation(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, eqn);
    //eqsys->SetInitialValue(OptionConstants::UQTY_E_FIELD, interp->Eval(t0), t0);
}

