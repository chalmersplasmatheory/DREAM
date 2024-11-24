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
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"

#include <algorithm>

using namespace DREAM;


#define MODULENAME "eqsys/j_ohm"


void SimulationGenerator::DefineOptions_j_ohm(Settings *s){
    s->DefineSetting(MODULENAME "/correctedConductivity", "Determines whether to use f_hot's natural ohmic current or the corrected (~Spitzer) value", (int_t) OptionConstants::CORRECTED_CONDUCTIVITY_ENABLED);
    s->DefineSetting(MODULENAME "/conductivityMode", "Determines which formula to use for the conductivity", (int_t) OptionConstants::CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS);

	DefineDataRT(MODULENAME, s, "data");
	DefineDataR(MODULENAME, s, "init");

	s->DefineSetting(MODULENAME "/Ip0", "When specifying an initial j_tot profile, re-scale the profile to get this total plasma current", (real_t)0.0);
	s->DefineSetting(MODULENAME "/j_type", "Current density profile type to initialize from.", (int_t)0.0);
}

/**
 * Construct the equation for the ohmic current density, 'j_ohm',
 * which represents j_\Omega / (B/Bmin) which is constant on the flux surface.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_ohm(
    EquationSystem *eqsys, Settings *s 
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    const len_t id_j_ohm = eqsys->GetUnknownID(OptionConstants::UQTY_J_OHM);
    const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
        
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
    std::string desc = "j_ohm = ";

    bool hottailMode = (eqsys->HasHotTailGrid()) && (eqsys->GetHotTailGrid()->GetNp2(0)==1);
    /** 
     * If using collfreq_mode FULL and not using hot tail mode (Nxi=1), 
     * calculate ohmic current by integrating the distribution, 
     * possibly with conductivity correction
     */
    if(eqsys->HasHotTailGrid() && (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) && !hottailMode){
        len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
        len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);

        // add total current carried by f_hot
        FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
        Op3->AddTerm(new CurrentDensityFromDistributionFunction(
                fluidGrid, eqsys->GetHotTailGrid(), id_j_ohm, id_f_hot,eqsys->GetUnknownHandler()
        ) );
        eqsys->SetOperator(id_j_ohm, id_f_hot, Op3);
        
        // subtract hot current (add with a scaleFactor of -1.0)
        FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
        Op4->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        desc += "integral(v_par*f_hot) - j_hot"; 
        eqsys->SetOperator(id_j_ohm, id_j_hot, Op4);
        
        OptionConstants::corrected_conductivity corrCond = (enum OptionConstants::corrected_conductivity)s->GetInteger(MODULENAME "/correctedConductivity");
        if(corrCond == OptionConstants::CORRECTED_CONDUCTIVITY_ENABLED){
            FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
            // add full ohmic current
            Op2->AddTerm(new CurrentFromConductivityTerm(
                            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()
            ) );
            // remove predicted numerical ohmic current (add with a scaleFactor of -1.0)
            Op2->AddTerm(new PredictedOhmicCurrentFromDistributionTerm(
                            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler(), -1.0
            ) );
            eqsys->SetOperator(id_j_ohm,id_E_field,Op2);
            desc += " + E*(sigma-sigma_num) [corrected]";

            // Initialization
            eqsys->initializer->AddRule(
                    id_j_ohm,
                    EqsysInitializer::INITRULE_EVAL_EQUATION,
                    nullptr,
                    // Dependencies
                    id_E_field,
                    id_f_hot,
                    EqsysInitializer::RUNAWAY_FLUID
            );

        } else {
            // Initialization
            eqsys->initializer->AddRule(
                    id_j_ohm,
                    EqsysInitializer::INITRULE_EVAL_EQUATION,
                    nullptr,
                    // Dependencies
                    id_f_hot
            );

        }
    // In all other cases, take full spitzer conductivity
    } else { 
        FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
        Op2->AddTerm(new CurrentFromConductivityTerm(
            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()
        ) );
        eqsys->SetOperator(id_j_ohm,id_E_field,Op2);
        desc += "sigma*E"; 
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
    eqsys->SetOperator(id_j_ohm,id_j_ohm,Op1,desc);    
}

