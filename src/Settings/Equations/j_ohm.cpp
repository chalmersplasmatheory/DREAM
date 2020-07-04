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
    s->DefineSetting(MODULENAME "/correctedConductivity", "Determines whether to use f_hot's natural ohmic current or the corrected (~Spitzer) value", (bool) true);
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
    const len_t id_j_ohm = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_J_OHM);
    const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");

    
    bool includeJHotInJOhm = ( eqsys->HasHotTailGrid() && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL));
    // If using full hot-tail grid, set ohmic current to the fast current
    if(includeJHotInJOhm){
        const len_t id_j_hot = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_J_HOT);
        FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
        FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
        Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        Op2->AddTerm(new FVM::IdentityTerm(fluidGrid));
        eqsys->SetOperator(id_j_ohm,id_j_ohm,Op1,"j_ohm = j_hot");
        eqsys->SetOperator(id_j_ohm,id_j_hot,Op2);

        // If also using the correctedConductivity setting, add the missing conductive current to j_ohm
        bool useCorrectedConductivity = (bool)s->GetBool(MODULENAME "/correctedConductivity");
        if(useCorrectedConductivity){
            FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
            
            // + (sigma-sigma_predicted)*E
            Op3->AddTerm(new CurrentFromConductivityTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()));
            Op3->AddTerm(new PredictedOhmicCurrentFromDistributionTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler(),-1.0));
            eqsys->SetOperator(id_j_ohm, id_E_field, Op3, "j_ohm = j_hot + Delta(sigma)*E");

            eqsys->initializer->AddRule(
                id_j_ohm,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                // Dependencies
                id_j_hot,
                id_E_field,
                EqsysInitializer::RUNAWAY_FLUID
            );

        } else {
            eqsys->initializer->AddRule(
                id_j_ohm,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                // Dependencies
                id_j_hot
            );
        }
    // Otherwise, calculate j_ohm from the full Sauter conductivity formula
    } else {
        FVM::Operator *eqn1 = new FVM::Operator(fluidGrid);
        FVM::Operator *eqn2 = new FVM::Operator(fluidGrid);
        
        // sigma*E
        eqn2->AddTerm(new CurrentFromConductivityTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()));

        // -j_ohm
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        
        eqsys->SetOperator(id_j_ohm, id_j_ohm, eqn1, "j_ohm = sigma*E");
        eqsys->SetOperator(id_j_ohm, id_E_field, eqn2);
        // Initialization
        eqsys->initializer->AddRule(
                id_j_ohm,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                // Dependencies
                id_E_field,
                EqsysInitializer::RUNAWAY_FLUID
        );
    }
}

