/**
 * Definition of equations relating to j_hot (the radial profile 
 * of parallel current density j_|| / (B/B_min) of hot electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "DREAM/Equations/Fluid/HotTailCurrentDensityFromDistributionFunction.hpp"
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
    EquationSystem *eqsys, Settings* s
) {
//    const real_t t0 = 0;

    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);
    len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    // If the hot-tail grid is enabled, we calculate j_hot as a
    // moment of the hot electron distribution function...
    if (hottailGrid) {
        if(hottailGrid->GetNp2(0)==1)
            ConstructEquation_j_hot_hottailMode(eqsys,s);
        else {
            len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

            FVM::Operator *eqn = new FVM::Operator(fluidGrid);

            CurrentDensityFromDistributionFunction *mq  = new CurrentDensityFromDistributionFunction(
                fluidGrid, hottailGrid, id_j_hot, id_f_hot
            );
            eqn->AddTerm(mq);
            eqsys->SetOperator(id_j_hot, id_f_hot, eqn, "Moment of f_hot - sigma_num*E");

            // Subtract predicted ohmic current (equal contribution is added to j_ohm)
            /*
            FVM::Operator *eqnE = new FVM::Operator(fluidGrid);
            eqnE->AddTerm(new PredictedOhmicCurrentFromDistributionTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler(),-1.0));
            eqsys->SetOperator(id_j_hot, id_E_field, eqnE);
            */

            // Identity part
            FVM::Operator *eqnIdent = new FVM::Operator(fluidGrid);
            eqnIdent->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
            eqsys->SetOperator(id_j_hot, id_j_hot, eqnIdent);

            // Initialize to zero
            //eqsys->SetInitialValue(OptionConstants::UQTY_N_HOT, nullptr, t0);

            // Set initialization method
            eqsys->initializer->AddRule(
                id_j_hot,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                // Dependencies
                id_f_hot,
                EqsysInitializer::RUNAWAY_FLUID,
                id_E_field
            );
        }

    // Otherwise, we set it to zero...
    } else {
        FVM::Operator *eqn = new FVM::Operator(fluidGrid);

        eqn->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
//        eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

        eqsys->SetOperator(id_j_hot, id_j_hot, eqn, "zero");
        // Set initialization method
        eqsys->initializer->AddRule(
            id_j_hot,
            EqsysInitializer::INITRULE_EVAL_EQUATION
        );
    }

}




/**
 * Sets j_hot = int( -E/nu_D * df_hot/dp, 0, p_cutoff ) + int( v*f_hot, p_cutoff, inf )
 * introducing a cutoff p_cutoff that ensures that the integrand is continuous.
 */
void SimulationGenerator::ConstructEquation_j_hot_hottailMode(
    EquationSystem *eqsys, Settings* /*s*/
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
//    eqsys->SetUnknown(OptionConstants::UQTY_J_HOT_P_CUT, fluidGrid);

    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    len_t id_jhot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_HOT);
    len_t id_fhot  = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
    len_t id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    len_t id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);


    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
    Op2->AddTerm(
        new HotTailCurrentDensityFromDistributionFunction(
            fluidGrid, eqsys->GetHotTailGrid(), unknowns,
            eqsys->GetHotTailCollisionHandler()->GetNuD()
        ) 
    );

    eqsys->SetOperator(id_jhot, id_jhot, Op1, "j_hot = int(df/dp) + int(f) [hot-tail mode]");
    eqsys->SetOperator(id_jhot, id_fhot, Op2);
    // Set initialization method
    eqsys->initializer->AddRule(
        id_jhot,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_fhot,
        id_Eterm,
        id_Tcold,
        EqsysInitializer::COLLQTYHDL_HOTTAIL
    );
}