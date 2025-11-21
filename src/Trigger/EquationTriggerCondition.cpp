
#include <DREAM/Trigger/EquationTriggerCondition.hpp>


using namespace DREAM;


/**
 * Constructor.
 */
EquationTriggerCondition::EquationTriggerCondition(
	FVM::Grid *g, FVM::UnknownQuantityHandler *u,
	const len_t nMultiples, bool saveTriggerState
) : grid(g), unknowns(u), nMultiples(nMultiples),
    saveTriggerState(saveTriggerState) {
	
	this->triggered = new bool[this->GetNCells()];
}

/**
 * Destructor.
 */
EquationTriggerCondition::~EquationTriggerCondition() {
	delete [] this->triggered;

	if (this->saveTriggerState)
		for (len_t i = 0; i < this->trigger_states.size(); i++)
			delete [] this->trigger_states[i];
}


/**
 * Check the trigger condition.
 */
void EquationTriggerCondition::CheckCondition(
	const real_t t, FVM::UnknownQuantityHandler *uqty
) {
	this->_CheckCondition(t, uqty);
	
	if (this->saveTriggerState)
		this->SaveTriggerState();
}


/**
 * Save the trigger state.
 */
void EquationTriggerCondition::SaveTriggerState() {
	this->trigger_states.push_back(this->triggered);
	this->triggered = new bool[this->GetNCells()];
}


