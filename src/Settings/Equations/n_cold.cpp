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
 * Define the equation for the cold electron density, 'n_cold'.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing hwo to construct
 *        the equation.
 */
void SimulationGenerator::DefineEquation_n_cold(
    EquationSystem *eqsys, Settings *s
) {
    enum uqty_n_cold_eqn eqn = (enum uqty_n_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (eqn) {
        case UQTY_N_COLD_EQN_PRESCRIBED:
            DefineEquation_n_cold_prescribed(eqsys, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for 'n_cold': %d.",
                eqn
            );
    }
}

/**
 * Define the equation describing a prescribed cold
 * electron density.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing hwo to construct
 *        the equation.
 */
void SimulationGenerator::DefineEquation_n_cold_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Equation *eqn = new FVM::Equation(eqsys->GetFluidGrid());

    enum FVM::Interpolator1D::interp_method im = FVM::Interpolator1D::INTERP_LINEAR;

    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), im);
    eqn->AddTerm(pp);

    //eqsys->SetEquation(UQTY_N_COLD, UQTY_N_COLD, eqn);
}

