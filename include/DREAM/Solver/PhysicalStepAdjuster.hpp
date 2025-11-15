#ifndef _DREAM_SOLVER_NONLINEAR_PHYSICAL_STEP_ADJUSTER_HPP
#define _DREAM_SOLVER_NONLINEAR_PHYSICAL_STEP_ADJUSTER_HPP

#include "DREAM/Solver/NewtonStepAdjuster.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
	class PhysicalStepAdjuster : public NewtonStepAdjuster {
	protected:
		IonHandler *ionHandler;

		// ID of unknown quantity which is currently
		// limited the step length the most.
		len_t limitingUnknown = 0;

		real_t MaximalStepLengthAtGridPoint(real_t, real_t, real_t);
		real_t MaximalPhysicalStepLength(const real_t*, const real_t*, len_t);

	public:
		PhysicalStepAdjuster(
			std::vector<len_t>& nu, FVM::UnknownQuantityHandler *uqh,
			IonHandler *ions
		);

		virtual real_t Adjust(
			len_t, const real_t*, const real_t*,
			Vec&, FVM::BlockMatrix*
		);
	};
}

#endif/*_DREAM_SOLVER_NONLINEAR_PHYSICAL_STEP_ADJUSTER_HPP*/
