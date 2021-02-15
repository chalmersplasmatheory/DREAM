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
#include <string>


using namespace DREAM;

#define MODULENAME "eqsys/n_hot"


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
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);

    // If the hot-tail grid is enabled, we calculate n_hot as a
    // moment of the hot electron distribution function...
    if (hottailGrid) {
        len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
        FVM::Operator *eqn = new FVM::Operator(fluidGrid);

        std::string desc = "integral(f_hot)";

        FVM::MomentQuantity::pThresholdMode pMode = (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
        real_t pThreshold = 0.0;
        enum OptionConstants::collqty_collfreq_mode collfreq_mode =
            (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
        if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
            // With collfreq_mode FULL, n_hot is defined as density above some threshold. 
            // For now: default definition of n_hot is p > 20*p_thermal 
            pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
            // TODO: format desc so that pThreshold is given explicitly (ie pThreshold*p_Te in this case) 
            desc = "integral(f_hot, p>pThreshold)"; 
        }
        eqn->AddTerm(new DensityFromDistributionFunction(
                fluidGrid, hottailGrid, id_n_hot, id_f_hot,eqsys->GetUnknownHandler(),pThreshold, pMode
            ));
        eqsys->SetOperator(id_n_hot, id_f_hot, eqn, desc);

        // Identity part
        FVM::Operator *eqnIdent = new FVM::Operator(fluidGrid);
        eqnIdent->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
        eqsys->SetOperator(id_n_hot, id_n_hot, eqnIdent);

        // Initialize from equation
        eqsys->initializer->AddRule(
            id_n_hot,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            // Dependencies
            id_f_hot,
            eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD)
        );

    // Otherwise, we set it to zero...
    } else {
        FVM::Operator *eqn = new FVM::Operator(fluidGrid);
        eqn->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqsys->SetOperator(id_n_hot, id_n_hot, eqn, "zero");

        // Initialization
        eqsys->initializer->AddRule(
            id_n_hot,
            EqsysInitializer::INITRULE_EVAL_EQUATION
        );
    }
}

