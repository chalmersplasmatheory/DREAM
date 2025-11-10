
#include <DREAM/Trigger/EquationTriggerCondition.hpp>


using namespace DREAM;


/**
 * Constructor.
 */
EquationTriggerCondition::EquationTriggerCondition(
	FVM::Grid *g, FVM::UnknownQuantityHandler *u,
	const len_t nMultiples
) : grid(g), unknowns(u), nMultiples(nMultiples) {
	
	this->triggered = new bool[this->GetNCells()];
}

/**
 * Destructor.
 */
EquationTriggerCondition::~EquationTriggerCondition() {
	delete [] this->triggered;
}


