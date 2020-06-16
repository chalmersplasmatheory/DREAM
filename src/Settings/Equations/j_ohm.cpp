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
#include "DREAM/Equations/Fluid/PredictedOhmicCurrentFromDistributionTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"

#include <algorithm>

using namespace DREAM;


#define MODULENAME "eqsys/j_ohm"


void SimulationGenerator::DefineOptions_j_ohm(Settings *s){
    s->DefineSetting(MODULENAME "/correctedConductivity", "Determines whether to use f_hot's natural ohmic current or the corrected (~Spitzer) value", (bool) false);
}

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

    
    bool useCorrectedConductivity = (bool)s->GetBool(MODULENAME "/correctedConductivity");
    if ((eqsys->HasHotTailGrid()) && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) && !useCorrectedConductivity ){
    /** 
     * Add the predicted ohmic current here, which was subtracted from j_hot
     */

        FVM::Operator *eqn1 = new FVM::Operator(fluidGrid);
        FVM::Operator *eqn2 = new FVM::Operator(fluidGrid);
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid, -1));
        eqn2->AddTerm(new PredictedOhmicCurrentFromDistributionTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()));
        eqsys->SetOperator(OptionConstants::UQTY_J_OHM,OptionConstants::UQTY_J_OHM,eqn1, "sigma_num*E");
        eqsys->SetOperator(OptionConstants::UQTY_J_OHM,OptionConstants::UQTY_E_FIELD,eqn2);

    // Otherwise, calculate it from the full Sauter conductivity formula
    } else {
        FVM::Operator *eqn1 = new FVM::Operator(fluidGrid);
        FVM::Operator *eqn2 = new FVM::Operator(fluidGrid);
        
        // sigma*E
        eqn2->AddTerm(new CurrentFromConductivityTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()));
        // -j_ohm
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        
        eqsys->SetOperator(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn1, "sigma_braams*Eterm");
        eqsys->SetOperator(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_E_FIELD, eqn2);


        
    }
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

