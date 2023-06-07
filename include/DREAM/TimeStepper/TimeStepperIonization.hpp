#ifndef _DREAM_TIME_STEPPER_IONIZATION_HPP
#define _DREAM_TIME_STEPPER_IONIZATION_HPP

#include <vector>
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
	class TimeStepperIonization : public TimeStepper {
	private:
		real_t tMax;
		real_t dt0, dtMax;
		real_t minSaveDt = 0, lastSaveTime = 0;

		real_t automaticstep = 1e-12;
		real_t safetyfactor = 50;

		real_t currentTime=0;
		real_t tscale0 = 0, dt;
		len_t id_n_cold;
		len_t currentStep = 0;

		len_t nr;
		real_t *timescales, *ncold;
		
		const len_t PROGRESSBAR_LENGTH = 80;
	public:
		TimeStepperIonization(
			const real_t tMax, const real_t dt0, const real_t dtMax,
			FVM::UnknownQuantityHandler*, const real_t automaticstep,
			const real_t safetyfactor, const real_t minSaveDt
		);
		~TimeStepperIonization();

		virtual real_t CurrentTime() const { return this->currentTime; }
		virtual bool IsFinished() override { return (this->currentTime >= this->tMax); }
		virtual bool IsSaveStep() override;
		virtual real_t NextTime() override;
		virtual void ValidateStep() override;

		real_t GetIonizationTimeScale();

		virtual void PrintProgress() override;
	};
}

#endif/*_DREAM_TIME_STEPPER_IONIZATION_HPP*/
