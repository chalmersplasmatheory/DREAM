/**
 * Definition of equations relating to n_hot (the radial
 * density of hot electrons).
 */
 
#include <sstream>
#include <iomanip>

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.enum.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include <string>


using namespace DREAM;


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
	const len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);

	ConstructEquation_n_hot_inner("eqsys/f_hot", eqsys, s, true);

	enum OptionConstants::eqn_trigger_type switchtype =
		(enum OptionConstants::eqn_trigger_type)s->GetInteger("eqsys/f_hot/switch/condition");
	
	// Set alternative equation?
	if (switchtype != OptionConstants::EQN_TRIGGER_TYPE_NONE) {
		eqsys->SetAssignToAlternativeEquation(id_n_hot, true);

		ConstructEquation_n_hot_inner("eqsys/f_hot/switch/equation", eqsys, s, false);

		EquationTriggerCondition *trig = LoadTriggerCondition(
			s, "eqsys/f_hot/switch", eqsys->GetFluidGrid(),
			eqsys->GetUnknownHandler()
		);
		eqsys->SetTriggerCondition(id_n_hot, trig);

		eqsys->SetAssignToAlternativeEquation(id_n_hot, false);
	}
}

void SimulationGenerator::ConstructEquation_n_hot_inner(
	const std::string &f_hot_modname, EquationSystem *eqsys,
	Settings *s, bool isMain
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);

    // If the hot-tail grid is enabled, we calculate n_hot as a
    // moment of the hot electron distribution function...
    if (hottailGrid) {
		enum OptionConstants::uqty_distribution_mode mode =
			(enum OptionConstants::uqty_distribution_mode)
				s->GetInteger(f_hot_modname + "/mode");

		if (mode == OptionConstants::UQTY_DISTRIBUTION_MODE_MAXWELLIAN) {
			len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
			FVM::Operator *eqn = new FVM::Operator(fluidGrid);

			eqn->AddTerm(new FVM::TransientTerm(fluidGrid, id_n_hot));
			eqsys->SetOperator(id_n_hot, id_n_hot, eqn, "dn_hot/dt = 0");

			real_t *n0 = _get_f_hot_data_r(s, "n0", fluidGrid->GetRadialGrid());
			eqsys->SetInitialValue(id_n_hot, n0);
			delete [] n0;
		} else {
			len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
			FVM::Operator *eqn = new FVM::Operator(fluidGrid);

			std::string desc = "integral(f_hot)";

			FVM::MomentQuantity::pThresholdMode pMode = (FVM::MomentQuantity::pThresholdMode)s->GetInteger(f_hot_modname + "/pThresholdMode");
			real_t pThreshold = 0.0;
			enum OptionConstants::collqty_collfreq_mode collfreq_mode =
				(enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
			if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
				// With collfreq_mode FULL, n_hot is defined as density above some threshold. 
				// For now: default definition of n_hot is p > 20*p_thermal 
				pThreshold = s->GetReal(f_hot_modname + "/pThreshold");
				
				std::ostringstream str;
				str <<std::fixed << std::setprecision(3) << pThreshold;
				switch(pMode){
					case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_MC:{
						desc = "integral(f_hot, p>"+str.str()+"*me*c)";
						break;
					} 
					case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL:{
						desc = "integral(f_hot, p>"+str.str()+"*pThermal)";
						break;
					}
					case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH:{
						desc = "integral(f_hot), smooth lower limit at p="+str.str()+"*pThermal";
						break;
					}
					case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_MC:{
						desc = "integral(f_hot, p<"+str.str()+"*me*c)";
						break;
					} 
					case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL:{
						desc = "integral(f_hot, p<"+str.str()+"*pThermal)";
						break;
					}
					case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH:{
						desc = "integral(f_hot), smooth upper limit at p="+str.str()+"*pThermal";
						break;
					}    
				}
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
			if (isMain) {
				eqsys->initializer->AddRule(
					id_n_hot,
					EqsysInitializer::INITRULE_EVAL_EQUATION,
					nullptr,
					// Dependencies
					id_f_hot,
					eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD)
				);
			}
		}
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

