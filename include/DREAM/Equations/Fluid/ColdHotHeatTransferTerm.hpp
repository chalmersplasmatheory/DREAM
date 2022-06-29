#ifndef _DREAM_EQUATIONS_FLUID_COLD_HOT_HEAT_TRANSFER_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_COLD_HOT_HEAT_TRANSFER_TERM_HPP

#include "FVM/Equation/MomentQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
	class ColdHotHeatTransferTerm : public FVM::MomentQuantity {
	private:
		real_t scaleFactor = 1.0;
		real_t *prevThresholdEnvelope = nullptr;

	public:
		ColdHotHeatTransferTerm(
			FVM::Grid*, FVM::Grid*, len_t, len_t,
			FVM::UnknownQuantityHandler*, real_t,
			FVM::MomentQuantity::pThresholdMode, real_t
		);
		virtual ~ColdHotHeatTransferTerm();

		virtual bool GridRebuilt() override;
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

		void StoreThresholdEnvelope();
	};
}

#endif/*_DREAM_EQUATIONS_FLUID_COLD_HOT_HEAT_TRANSFER_TERM_HPP*/
