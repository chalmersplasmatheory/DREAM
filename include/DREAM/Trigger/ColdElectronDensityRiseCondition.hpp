#ifndef _DREAM_TRIGGER_COLD_ELECTRON_DENSITY_RISE_CONDITION_HPP
#define _DREAM_TRIGGER_COLD_ELECTRON_DENSITY_RISE_CONDITION_HPP

#include <DREAM/Trigger/EquationTriggerCondition.hpp>

namespace DREAM {
	class ColdElectronDensityRiseCondition : public EquationTriggerCondition {
	protected:
		real_t sensitivity = 0.01;

		len_t id_n_cold, id_n_hot;

	public:
		ColdElectronDensityRiseCondition(
			FVM::Grid*, FVM::UnknownQuantityHandler*, const real_t
		);
		virtual ~ColdElectronDensityRiseCondition();

		virtual void CheckCondition(FVM::UnknownQuantityHandler*) override;
	};
}

#endif/*_DREAM_TRIGGER_COLD_ELECTRON_DENSITY_RISE_CONDITION_HPP*/
