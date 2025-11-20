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
		real_t maximalPhysicalStepLength = 1;
		bool hasAdjusted = false;

		std::vector<len_t> ids_nonNegativeQuantities;
		len_t id_ni;

		real_t MaximalStepLengthAtGridPoint(real_t, real_t, real_t);
		real_t MaximalPhysicalStepLength(const real_t*, const real_t*, len_t);

	public:
		PhysicalStepAdjuster(
			std::vector<len_t>& nu, FVM::UnknownQuantityHandler *uqh,
			IonHandler *ions, const len_t
		);
		virtual ~PhysicalStepAdjuster();

		virtual bool AdjustmentNeeded(const len_t, Vec&, FVM::BlockMatrix*) override;
		virtual void AdjustSolution(const len_t, real_t*) override;
		virtual void Reset(Vec&, FVM::BlockMatrix*) override;

		virtual real_t GetCurrentDamping() override { return this->maximalPhysicalStepLength; }

		virtual void SetX0(const len_t, const real_t*, const real_t*) override;
	};
}

#endif/*_DREAM_SOLVER_NONLINEAR_PHYSICAL_STEP_ADJUSTER_HPP*/
