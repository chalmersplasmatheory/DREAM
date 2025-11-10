#ifndef _DREAM_TRIGGER_EQUATION_TRIGGER_CONDITION_HPP
#define _DREAM_TRIGGER_EQUATION_TRIGGER_CONDITION_HPP

#include <FVM/UnknownQuantityHandler.hpp>

namespace DREAM {
	class EquationTriggerCondition {
	protected:
		FVM::Grid *grid;
		FVM::UnknownQuantityHandler *unknowns;
		const len_t nMultiples = 1;

		bool *triggered;

		void SetTriggered(const len_t i, bool v) { this->triggered[i] = v; }

	public:
		EquationTriggerCondition(FVM::Grid*, FVM::UnknownQuantityHandler*, const len_t nMultiples);
		virtual ~EquationTriggerCondition();

		virtual void CheckCondition(const real_t, FVM::UnknownQuantityHandler*) = 0;
		const bool *GetTriggerMask() { return this->triggered; }
		const len_t GetNCells() { return this->nMultiples * this->grid->GetNCells(); }

		bool IsTriggered(const len_t i) { return this->triggered[i]; }
	};
}

#endif/*_DREAM_TRIGGER_EQUATION_TRIGGER_CONDITION_HPP*/
