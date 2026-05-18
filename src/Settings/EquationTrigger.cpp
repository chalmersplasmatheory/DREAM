/**
 * Helper routines for trigger conditions.
 */

#include <string>
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

#include "DREAM/Trigger/ColdElectronDensityRiseCondition.hpp"
#include "DREAM/Trigger/EquationTriggerCondition.hpp"
#include "DREAM/Trigger/TimeTrigger.hpp"


using namespace DREAM;
using namespace std;


/**
 * Define general options for a trigger condition.
 */
void SimulationGenerator::DefineOptions_TriggerCondition(
	Settings *s, const string& modulename
) {
	s->DefineSetting(modulename + "/condition", "Type of condition to use for triggering the equation switch", (int_t)OptionConstants::EQN_TRIGGER_TYPE_NONE);
	s->DefineSetting(modulename + "/trigger_time", "For the 'TIME' equation trigger: simulation time at which to trigger the condition", (real_t)0);
	s->DefineSetting(modulename + "/sensitivity", "For the 'COLD_ELECTRON_RISE' equation trigger: ratio of n_cold/n_hot at which to trigger the alternative equation", (real_t)0.01);
}

/**
 * Load the trigger condition for the specified module.
 */
EquationTriggerCondition *SimulationGenerator::LoadTriggerCondition(
	Settings *s, const string& name, FVM::Grid *grid,
	FVM::UnknownQuantityHandler *uqh, const len_t nMultiples
) {
	enum OptionConstants::eqn_trigger_type ttype =
		(enum OptionConstants::eqn_trigger_type)s->GetInteger(name + "/condition");
	
	EquationTriggerCondition *cond;
	switch (ttype) {
		case OptionConstants::EQN_TRIGGER_TYPE_NONE:
			cond = nullptr;
			break;

		case OptionConstants::EQN_TRIGGER_TYPE_TIME: {
			real_t ttime = s->GetReal(name + "/trigger_time");
			cond = new TimeTrigger(grid, uqh, nMultiples, ttime);
		} break;

		case OptionConstants::EQN_TRIGGER_TYPE_COLD_ELECTRON_DENSITY_RISE: {
			real_t sensitivity = s->GetReal(name + "/sensitivity");
			cond = new ColdElectronDensityRiseCondition(grid, uqh, nMultiples, sensitivity);
		} break;
		
		default:
			throw DREAMException("Unrecognized equation trigger condition: %d.", ttype);
	}

	return cond;
}


