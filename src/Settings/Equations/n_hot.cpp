/**
 * Definition of equations relating to n_hot (the radial
 * density of hot electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "equationsystem/n_hot"


/**
 * Construct the equation for the hot electron density, 'n_hot'.
 * If the hot-tail grid is enabled, n_hot will be the integral of
 * the hot electron distribution. If it does not exist, we set the
 * number of hot electrons to zero and call electrons leaving the
 * grid "runaways" instead (and, thus, counting them with 'n_re').
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_n_hot(
    EquationSystem *eqsys, Settings*
) {
    const real_t t0 = 0;

    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
    len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

    // If the hot-tail grid is enabled, we calculate n_hot as a
    // moment of the hot electron distribution function...
    if (hottailGrid) {
        FVM::Equation *eqn = new FVM::Equation(fluidGrid);

        DensityFromDistributionFunction *mq  = new DensityFromDistributionFunction(
            fluidGrid, hottailGrid, id_n_hot, id_f_hot
        );
        eqn->AddTerm(mq);
        eqsys->SetEquation(id_n_hot, id_f_hot, eqn);

        // Identity part
        FVM::Equation *eqnIdent = new FVM::Equation(fluidGrid);
        eqnIdent->AddTerm(new FVM::IdentityTerm(fluidGrid));
        eqsys->SetEquation(id_n_hot, id_n_hot, eqnIdent);

        // Initialize to zero
        //eqsys->SetInitialValue(OptionConstants::UQTY_N_HOT, nullptr, t0);
    // Otherwise, we set it to zero...
    } else {
        FVM::Equation *eqn = new FVM::Equation(fluidGrid);

        FVM::ConstantParameter *np = new FVM::ConstantParameter(fluidGrid, 0);
        eqn->AddTerm(np);

        eqsys->SetEquation(id_n_hot, id_n_hot, eqn);
    }
}

