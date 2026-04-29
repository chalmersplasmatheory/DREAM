/**
 * A simple equation switch trigger condition which triggers after a given time.
 */

#include "DREAM/Trigger/TimeTrigger.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
TimeTrigger::TimeTrigger(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *uqh,
	const len_t nMultiples, const real_t trigger_time,
	bool saveTriggerState
) : EquationTriggerCondition(grid, uqh, nMultiples, saveTriggerState), trigger_time(trigger_time) { }


/**
 * Destructor.
 */
TimeTrigger::~TimeTrigger() {}


/**
 * Check if the condition is triggered.
 */
void TimeTrigger::CheckCondition(const real_t t, FVM::UnknownQuantityHandler*) {
	const len_t N = this->GetNCells();
	bool v = (t >= this->trigger_time);
	for (len_t i = 0; i < N; i++)
		triggered[i] = v;
}


