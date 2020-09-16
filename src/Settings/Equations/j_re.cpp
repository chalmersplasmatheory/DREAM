/**
 * Definition of equations relating to j_hot (the radial profile 
 * of parallel current density j_|| / (B/B_min) of hot electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/j_re"


/**
 * Construct the equation for the hot parallel current, 'j_hot'.
 * If the hot-tail grid is enabled, j_hot will be an integral of
 * the hot electron distribution. If it does not exist, it is set
 * to 0.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_re(
    EquationSystem *eqsys, Settings* /*s*/
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    len_t id_j_re = eqsys->GetUnknownID(OptionConstants::UQTY_J_RE);
    len_t id_n_re = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);


    // Identity part
    FVM::Operator *eqnIdent = new FVM::Operator(fluidGrid);
    eqnIdent->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
    eqsys->SetOperator(id_j_re, id_j_re, eqnIdent);

    FVM::Operator *eqn = new FVM::Operator(fluidGrid);

    // if runawayGrid is enabled, take moment of f_re, otherwise e*c*n_re
    if (runawayGrid) {
        len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
        eqn->AddTerm(new CurrentDensityFromDistributionFunction(
            fluidGrid, runawayGrid, id_j_re, id_f_re
        ));
        eqsys->SetOperator(id_j_re, id_f_re, eqn, "Moment of f_re");

        // Set initialization method
        eqsys->initializer->AddRule(
            id_j_re,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            // Dependencies
            id_f_re
        );

    // Otherwise, we set it to zero...
    } else {
        eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, Constants::ec * Constants::c));
        eqsys->SetOperator(id_j_re, id_n_re, eqn, "j_re = e*c*n_re");
        // Set initialization method
        eqsys->initializer->AddRule(
            id_j_re,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            id_n_re
        );
    }

}
