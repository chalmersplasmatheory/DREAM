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
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);
    len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    // If the hot-tail grid is enabled, we calculate j_hot as a
    // moment of the hot electron distribution function...
    if (hottailGrid) {
        if(hottailGrid->GetNp2(0)==1) // XXX: assumes we don't switch mode between radii
            ConstructEquation_j_hot_hottailMode(eqsys,s);
        else {
            len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

            // Identity part
            FVM::Operator *Op0 = new FVM::Operator(fluidGrid);
            Op0->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
            eqsys->SetOperator(id_j_hot, id_j_hot, Op0);

            std::string desc = "Moment of f_hot";

            FVM::MomentQuantity::pThresholdMode pMode = (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
            real_t pThreshold = 0.0;
            enum OptionConstants::collqty_collfreq_mode collfreq_mode =
                (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
            if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
                // With collfreq_mode FULL, n_hot is defined as density above some threshold.
                pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
                // TODO: format desc so that pThreshold is given explicitly (ie 5*p_Te in this case) 
                desc = "integral(v_par*f_hot, p>pThreshold)"; 
            }
            FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
            Op1->AddTerm(new CurrentDensityFromDistributionFunction(
                    fluidGrid, hottailGrid, id_j_hot, id_f_hot,eqsys->GetUnknownHandler(),pThreshold, pMode)
                );
            eqsys->SetOperator(id_j_hot, id_f_hot, Op1, desc);

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
        eqsys->SetOperator(id_j_hot, id_j_hot, eqn, "zero");
        // Set initialization method
        eqsys->initializer->AddRule(
            id_j_hot,
            EqsysInitializer::INITRULE_EVAL_EQUATION
        );
    }

}


/**
 * Based on the two formulas
 *      j1 = int( -E/nu_D * df_hot/dp , dp ) [Lorentz limit]
 *      j2 = int( v*f_hot * sign(E) , dp )   [theta << 1 limit],
 * constructs
 *      j_hot = int( dj1*dj2/sqrt(dj1^2+dj2^2), dp)
 * as roughly the smallest of these (Lorentz limit for 
 * low p or weak E, otherwise <v_par> ~ <v>).
 */
void SimulationGenerator::ConstructEquation_j_hot_hottailMode(
    EquationSystem *eqsys, Settings* s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    len_t id_jhot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_HOT);
    len_t id_fhot  = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
    len_t id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    len_t id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

    enum OptionConstants::collqty_collfreq_mode collfreq_mode = (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    Op2->AddTerm(
        new HotTailCurrentDensityFromDistributionFunction(
            fluidGrid, eqsys->GetHotTailGrid(), unknowns,
            eqsys->GetHotTailCollisionHandler()->GetNuD(),
            collfreq_mode
        ) 
    );

    eqsys->SetOperator(id_jhot, id_jhot, Op1, "j_hot = int(-E/nu_D * df_hot/dp) + int(v*f_hot) [hot-tail mode]");
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