/**
 * Definition of equations relating to j_hot (the radial profile 
 * of parallel current density j_|| / (B/B_min) of hot electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "DREAM/Equations/Fluid/PredictedOhmicCurrentFromDistributionTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/j_hot"


/**
 * Construct the equation for the hot parallel current, 'j_hot'.
 * If the hot-tail grid is enabled, j_hot will be an integral of
 * the hot electron distribution. If it does not exist, it is set
 * to 0.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_hot(
    EquationSystem *eqsys, Settings*
) {
//    const real_t t0 = 0;

    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);
    len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
    len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    // If the hot-tail grid is enabled, we calculate j_hot as a
    // moment of the hot electron distribution function...
    if (hottailGrid) {
        FVM::Operator *eqn = new FVM::Operator(fluidGrid);

        CurrentDensityFromDistributionFunction *mq  = new CurrentDensityFromDistributionFunction(
            fluidGrid, hottailGrid, id_j_hot, id_f_hot
        );
        eqn->AddTerm(mq);
        eqsys->SetOperator(id_j_hot, id_f_hot, eqn, "Moment of f_hot - sigma_num*E");

        // Subtract predicted ohmic current (equal contribution is added to j_ohm)
        FVM::Operator *eqnE = new FVM::Operator(fluidGrid);
        eqnE->AddTerm(new PredictedOhmicCurrentFromDistributionTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler(),-1.0));
        eqsys->SetOperator(id_j_hot, id_E_field, eqnE);
        // Identity part
        FVM::Operator *eqnIdent = new FVM::Operator(fluidGrid);
        eqnIdent->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
        eqsys->SetOperator(id_j_hot, id_j_hot, eqnIdent);

        // Initialize to zero
        //eqsys->SetInitialValue(OptionConstants::UQTY_N_HOT, nullptr, t0);
    // Otherwise, we set it to zero...
    } else {
        FVM::Operator *eqn = new FVM::Operator(fluidGrid);

        eqn->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

        eqsys->SetOperator(id_j_hot, id_j_hot, eqn, "zero");
    }

    // Set initialization method
    const len_t id_n_cold  = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const len_t id_n_i     = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    const len_t id_T_cold  = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);
    eqsys->initializer->AddRule(
        id_j_hot,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_f_hot,id_E_field, id_n_cold, id_n_i, id_T_cold
    );

}

