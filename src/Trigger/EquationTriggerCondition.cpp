
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
 * Save the trigger state.
 */
void EquationTriggerCondition::SaveTriggerState() {
	if (!this->saveTriggerState)
		return;

	const len_t N = this->GetNCells();
	bool *s = new bool[N];
	for (len_t i = 0; i < N; i++)
		s[i] = this->triggered[i];

	this->trigger_states.push_back(s);
}


