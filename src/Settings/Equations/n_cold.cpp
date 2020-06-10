/**
 * Definition of equations relating to n_cold.
 */

#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/NColdFromQuasiNeutrality.hpp"
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
    FVM::Equation *eqn = new FVM::Equation(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetEquation(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_COLD, eqn, "Prescribed");
    //eqsys->SetInitialValue(OptionConstants::UQTY_N_COLD, interp->Eval(t0), t0);
}

/**
 * Construct the equation describing the cold electron density as
 * (if superthermal grid and hot electrons are transferred to n_cold via lower p boundary:)
 *   n_cold = n_free - n_hot - n_re
 * (if particle conserving collision operator:)
 *   n_cold = n_hot + n_re
 */
void SimulationGenerator::ConstructEquation_n_cold_selfconsistent(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    const len_t id_nhot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
    const len_t id_nre  = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);

    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");

    if (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL) {
        FVM::Equation *eqn = new FVM::Equation(fluidGrid);

        eqn->AddTerm(new NColdFromQuasiNeutrality(fluidGrid, eqsys->GetIonHandler(), id_nhot, id_nre));
        eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

        eqsys->SetEquation(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_COLD, eqn, "Self-consistent");
    } else {
        FVM::Equation *eqn0 = new FVM::Equation(fluidGrid);
        FVM::Equation *eqn1 = new FVM::Equation(fluidGrid);
        FVM::Equation *eqn2 = new FVM::Equation(fluidGrid);

        eqn0->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid));
        eqn2->AddTerm(new FVM::IdentityTerm(fluidGrid));

        eqsys->SetEquation(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_COLD, eqn0, "n_cold = n_hot + n_re");
        eqsys->SetEquation(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_HOT, eqn1);
        eqsys->SetEquation(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_RE, eqn2);
    }
}

