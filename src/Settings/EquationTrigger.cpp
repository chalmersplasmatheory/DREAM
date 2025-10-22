/**
 * Helper routines for trigger conditions.
 */

#include <string>
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

#include "DREAM/Trigger/ColdElectronDensityRiseCondition.hpp"
#include "DREAM/Trigger/EquationTriggerCondition.hpp"


using namespace DREAM;
using namespace std;


EquationTriggerCondition *SimulationGenerator::LoadTriggerCondition(
	Settings *s, const string& name, FVM::Grid *grid,
	FVM::UnknownQuantityHandler *uqh
) {
	enum OptionConstants::eqn_trigger_type ttype =
		(enum OptionConstants::eqn_trigger_type)s->GetInteger(name + "/condition");
	
	EquationTriggerCondition *cond;
	switch (ttype) {
		case OptionConstants::EQN_TRIGGER_TYPE_NONE:
			cond = nullptr;
			break;

		case OptionConstants::EQN_TRIGGER_TYPE_COLD_ELECTRON_DENSITY_RISE: {
			real_t sensitivity = s->GetReal(name + "/sensitivity");
			cond = new ColdElectronDensityRiseCondition(grid, uqh, sensitivity);
		} break;
		
		default:
			throw DREAMException("Unrecognized equation trigger condition: %d.", ttype);
	}

	return cond;
}


