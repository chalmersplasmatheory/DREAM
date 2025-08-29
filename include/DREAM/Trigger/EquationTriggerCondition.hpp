#ifndef _DREAM_TRIGGER_EQUATION_TRIGGER_CONDITION_HPP
#define _DREAM_TRIGGER_EQUATION_TRIGGER_CONDITION_HPP

#include <FVM/UnknownQuantityHandler.hpp>

namespace DREAM {
	class EquationTriggerCondition {
	protected:
		FVM::Grid *grid;
		FVM::UnknownQuantityHandler *unknowns;
		bool *triggered;

		void SetTriggered(bool v) { this->triggered = v; }

	public:
		EquationTriggerCondition(FVM::Grid*, FVM::UnknownQuantityHandler*);
		virtual ~EquationTriggerCondition();

		virtual void CheckCondition(UnknownQuantityHandler*) = 0;

		bool IsTriggered(const len_t i) { return this->triggered[i]; }
	};
}

#endif/*_DREAM_TRIGGER_EQUATION_TRIGGER_CONDITION_HPP*/
