/**
 * Definition of equations relating to the radial profile 
 * of ohmic current density j_ohm.
 * The quantity j_ohm corresponds to 
 * j_\Omega / (B/Bmin),
 * which is constant on flux surfaces and proportional to
 * j_ohm ~ \sigma*E_term,
 * where sigma is a neoclassical conductivity
 * containing various geometrical corrections
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentFromConductivityTerm.hpp"

#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"

#include <algorithm>

using namespace DREAM;


#define MODULENAME "eqsys/j_ohm"


/**
 * Construct the equation for the ohmic current density, 'j_ohm',
 * which represents j_\Omega / (B/Bmin) which is constant on the flux surface.
 * This is zero when the hot tail grid uses collfreqmode FULL,
 * in which case the ohmic current is part of f_hot. 
 * The ohmic current j_ohm is calculated from a semi-analytical 
 * electric conductivity formula.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_ohm(
    EquationSystem *eqsys, Settings *s 
) {

    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();

    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");

    
    if ((eqsys->HasHotTailGrid()) && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)) {
    /** 
     * TODO: If not full momentum-conserving operator, this would be a good place to 
     *       insert Linnea's "corrected conductivity" routine. For now setting to 0.
     */
        FVM::Equation *eqn = new FVM::Equation(fluidGrid);

        eqn->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        //eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn, "zero");

        // Initialization
        eqsys->initializer->AddRule(
            OptionConstants::UQTY_J_OHM,
            EqsysInitializer::INITRULE_EVAL_EQUATION
        );
    // Otherwise, calculate it from conductivity
    } else {
        FVM::Equation *eqn1 = new FVM::Equation(fluidGrid);
        FVM::Equation *eqn2 = new FVM::Equation(fluidGrid);
        
        // sigma*E
        eqn2->AddTerm(new CurrentFromConductivityTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()));
        // -j_ohm
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn1, "j_ohm = sigma*Eterm");
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_E_FIELD, eqn2);

        // Initialization
        const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);
        // (conductivity depends on these)
        const len_t id_n_cold  = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
        const len_t id_n_i     = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
        const len_t id_T_cold  = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);

        eqsys->initializer->AddRule(
            OptionConstants::UQTY_J_OHM,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            // Dependencies
            id_E_field, id_n_cold, id_n_i, id_T_cold
        );
    }
}

