/**
 * Definition of equations relating to j_hot (the radial profile 
 * of parallel current density j_|| / (B/B_min) of hot electrons).
 */
 
#include <sstream>
#include <iomanip>

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "DREAM/Equations/Fluid/CurrentFromHotConductivityTerm.hpp"
#include "DREAM/Equations/Fluid/HotTailCurrentDensityFromDistributionFunction.hpp"
#include "DREAM/Equations/Fluid/PredictedOhmicCurrentFromDistributionTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/j_hot"


/**
 * Define options for the hot electron current density module.
 */
void SimulationGenerator::DefineOptions_j_hot(Settings *s) {
	DefineOptions_j_hot_inner(s, MODULENAME);
	DefineOptions_TriggerCondition(s, MODULENAME "/switch");
	DefineOptions_j_hot_inner(s, MODULENAME "/switch/equation");
}

void SimulationGenerator::DefineOptions_j_hot_inner(Settings *s, const std::string &modulename) {
	s->DefineSetting(modulename+ "/type", "Type of equation to use for the hot electron current density.", (int_t)OptionConstants::UQTY_J_HOT_EQN_MOMENT);
}


/**
 * Construct the equation for the hot electron current density.
 */
void SimulationGenerator::ConstructEquation_j_hot(
    EquationSystem *eqsys, Settings *s
) {
	const len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);

	// Construct main equation
	ConstructEquation_j_hot_inner(MODULENAME, eqsys, s);

	enum OptionConstants::eqn_trigger_type switchtype =
		(enum OptionConstants::eqn_trigger_type)s->GetInteger(MODULENAME "/switch/condition");
	
	// Set alternative equation?
	if (switchtype != OptionConstants::EQN_TRIGGER_TYPE_NONE) {
		eqsys->SetAssignToAlternativeEquation(id_j_hot, true);

		ConstructEquation_j_hot_inner(MODULENAME "/switch/equation", eqsys, s);

		EquationTriggerCondition *trig = LoadTriggerCondition(s, MODULENAME "/switch", eqsys->GetFluidGrid(), eqsys->GetUnknownHandler());
		eqsys->SetTriggerCondition(id_j_hot, trig);

		eqsys->SetAssignToAlternativeEquation(id_j_hot, false);
	}
}

void SimulationGenerator::ConstructEquation_j_hot_inner(
	const std::string &modulename, EquationSystem *eqsys,
	Settings *s
) {
	enum OptionConstants::uqty_j_hot_eqn eqn_type =
		(enum OptionConstants::uqty_j_hot_eqn)s->GetInteger(modulename + "/type");
	
	switch (eqn_type) {
		case OptionConstants::UQTY_J_HOT_EQN_MOMENT:
			ConstructEquation_j_hot_moment(eqsys, s);
			break;

		case OptionConstants::UQTY_J_HOT_EQN_OHMIC:
			ConstructEquation_j_hot_ohmic(eqsys, s);
			break;

		default:
			throw SettingsException(
				"Unrecognized equation type for 'j_hot': %d.",
				eqn_type
			);
	}
}

/**
 * Construct the equation for the hot parallel current, 'j_hot'.
 * With this option, we take the hot current to be governed by
 * Ohm's law, evaluated at the hot electron temperature and
 * density.
 */
void SimulationGenerator::ConstructEquation_j_hot_ohmic(
	EquationSystem *eqsys, Settings*
) {
	FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();

	const len_t id_j_hot   = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);
	const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);

	FVM::Operator *Op_j = new FVM::Operator(fluidGrid);
	FVM::Operator *Op_E = new FVM::Operator(fluidGrid);

	Op_j->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
	Op_E->AddTerm(new CurrentFromHotConductivityTerm(
		fluidGrid, eqsys->GetUnknownHandler(),
		eqsys->GetREFluid(), eqsys->GetIonHandler()
	));

	eqsys->SetOperator(id_j_hot, id_j_hot, Op_j, "j_hot = sigma*E");
	eqsys->SetOperator(id_j_hot, id_E_field, Op_E);

	eqsys->initializer->AddRule(
		id_j_hot,
		EqsysInitializer::INITRULE_EVAL_EQUATION,
		nullptr,
		// Dependencies
		id_E_field,
		EqsysInitializer::RUNAWAY_FLUID
	);
}

/**
 * Construct the equation for the hot parallel current, 'j_hot'.
 * If the hot-tail grid is enabled, j_hot will be an integral of
 * the hot electron distribution. If it does not exist, it is set
 * to 0.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_hot_moment(
	EquationSystem *eqsys, Settings *s
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

            FVM::MomentQuantity::pThresholdMode pMode = (FVM::MomentQuantity::pThresholdMode)_get_f_hot_int(s, "pThresholdMode");
            real_t pThreshold = 0.0;
            enum OptionConstants::collqty_collfreq_mode collfreq_mode =
                (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
            if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
                // With collfreq_mode FULL, n_hot is defined as density above some threshold.
                pThreshold = (real_t)_get_f_hot_real(s, "pThreshold"); 
                
                std::ostringstream str;
                str <<std::fixed << std::setprecision(3) << pThreshold;
                switch(pMode){
                    case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_MC:{
                        desc = "integral(v_par*f_hot, p>"+str.str()+"*me*c)";
                        break;
                    } 
                    case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL:{
                        desc = "integral(v_par*f_hot, p>"+str.str()+"*pThermal)";
                        break;
                    }
                    case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH:{
                        desc = "integral(v_par*f_hot), smooth lower limit at p="+str.str()+"*pThermal";
                        break;
                    }
                    case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_MC:{
                        desc = "integral(v_par*f_hot, p<"+str.str()+"*me*c)";
                        break;
                    } 
                    case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL:{
                        desc = "integral(v_par*f_hot, p<"+str.str()+"*pThermal)";
                        break;
                    }
                    case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH:{
                        desc = "integral(v_par*f_hot), smooth upper limit at p="+str.str()+"*pThermal";
                        break;
                    }    
                }
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
    bool withFullIonJacobian = (bool)_get_f_hot_bool(s, "fullIonJacobian");
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    Op2->AddTerm(
        new HotTailCurrentDensityFromDistributionFunction(
            fluidGrid, eqsys->GetHotTailGrid(), unknowns,
            eqsys->GetHotTailCollisionHandler()->GetNuD(),
            collfreq_mode, withFullIonJacobian
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

/**
 * Returns the full path to the appropriate quantity 'name' to use.
 * This quantity will either be located under '/eqsys/f_hot' or
 * '/eqsys/f_hot/switch/equation', depending on whether it is set for the
 * main or alternative equations.
 */
std::string SimulationGenerator::_get_f_hot_subgroup(Settings *s, const std::string& name) {
	enum OptionConstants::eqn_trigger_type switchtype =
		(enum OptionConstants::eqn_trigger_type)s->GetInteger("eqsys/f_hot/switch/condition");
	
	// Main equation?
	if (switchtype == OptionConstants::EQN_TRIGGER_TYPE_NONE)
		return "eqsys/f_hot/" + name;
	else
		return "eqsys/f_hot/switch/equation/" + name;
}
bool SimulationGenerator::_get_f_hot_bool(Settings *s, const std::string& name) {
	const std::string &path = _get_f_hot_subgroup(s, name);
	return s->GetBool(path);
}
int_t SimulationGenerator::_get_f_hot_int(Settings *s, const std::string& name) {
	const std::string &path = _get_f_hot_subgroup(s, name);
	return s->GetInteger(path);
}
real_t SimulationGenerator::_get_f_hot_real(Settings *s, const std::string& name) {
	const std::string &path = _get_f_hot_subgroup(s, name);
	return s->GetReal(path);
}

