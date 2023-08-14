#ifndef _DREAM_FVM_EQUATION_PRESCRIBED_KINETIC_PARAMETER_HPP
#define _DREAM_FVM_EQUATION_PRESCRIBED_KINETIC_PARAMETER_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "DREAM/Settings/LoadData.hpp"

namespace DREAM {
	class PrescribedKineticParameter : public FVM::PrescribedParameter {
	private:
		struct dream_4d_data *data = nullptr;

		real_t *evalData = nullptr;
		len_t prevTimestep = 0;
	
	public:
		PrescribedKineticParameter(FVM::Grid*);
		PrescribedKineticParameter(FVM::Grid*, struct dream_4d_data*);
		~PrescribedKineticParameter();

		void DeallocateData();

		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
	};
}

#endif/*_DREAM_FVM_EQUATION_PRESCRIBED_KINETIC_PARAMETER_HPP*/

