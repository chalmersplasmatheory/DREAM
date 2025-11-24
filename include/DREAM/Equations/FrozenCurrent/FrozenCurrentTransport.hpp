#ifndef _DREAM_EQUATIONS_KINETIC_FROZEN_CURRENT_TRANSPORT_HPP
#define _DREAM_EQUATIONS_KINETIC_FROZEN_CURRENT_TRANSPORT_HPP

#include "FVM/Equation/DiffusionTerm.hpp"


namespace DREAM {
	class FrozenCurrentTransport : public FVM::DiffusionTerm {
	public:
		enum TransportMode {
			TRANSPORT_MODE_CONSTANT,	// Transport independent of momentum
			TRANSPORT_MODE_BETAPAR		// Transport scales as normalized parallel speed
		};
	protected:
		enum TransportMode transportMode;
		len_t id_D_I;

	public:
		FrozenCurrentTransport(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			enum TransportMode, bool allocCoefficients=false
		);
		virtual ~FrozenCurrentTransport();

		virtual bool GridRebuilt() override;
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

		virtual void SetPartialDiffusionTerm(len_t, len_t) override;
	};
}

#endif/*_DREAM_EQUATIONS_KINETIC_FROZEN_CURRENT_TRANSPORT_HPP*/
