
#include <DREAM/Trigger/EquationTriggerCondition.hpp>


using namespace DREAM;


/**
 * Constructor.
 */
EquationTriggerCondition::EquationTriggerCondition(
	FVM::Grid *g, FVM::UnknownQuantityHandler *u
) : grid(g), unknowns(u) {
	
	this->triggered = new bool[g->GetNCells()];
}

/**
 * Destructor.
 */
EquationTriggerCondition::~EquationTriggerCondition() {
	delete [] this->triggered;
}


