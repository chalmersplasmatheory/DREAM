#ifndef _DREAM_TRIGGER_TIME_TRIGGER_HPP
#define _DREAM_TRIGGER_TIME_TRIGGER_HPP

#include <DREAM/Trigger/EquationTriggerCondition.hpp>

namespace DREAM {
	class TimeTrigger : public EquationTriggerCondition {
	protected:
		real_t trigger_time = 0;

		len_t id_n_cold, id_n_hot;

	public:
		TimeTrigger(
			FVM::Grid*, FVM::UnknownQuantityHandler*, const len_t,
			const real_t, bool saveTriggerState=true
		);
		virtual ~TimeTrigger();

		virtual void CheckCondition(const real_t, FVM::UnknownQuantityHandler*) override;
	};
}

#endif/*_DREAM_TRIGGER_TIME_TRIGGER_HPP*/
