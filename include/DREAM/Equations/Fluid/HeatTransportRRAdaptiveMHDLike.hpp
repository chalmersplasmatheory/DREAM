#ifndef _DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RR_ADAPTIVE_MHD_LIKE_HPP
#define _DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RR_ADAPTIVE_MHD_LIKE_HPP

#include "DREAM/Equations/AdaptiveMHDLikeTransportTerm.hpp"
#include "DREAM/Equations/Fluid/HeatTransportRechesterRosenbluth.hpp"


namespace DREAM {
	class HeatTransportRRAdaptiveMHDLike :
		public AdaptiveMHDLikeTransportTerm,
		public HeatTransportRechesterRosenbluth {
	protected:
		real_t dBOverB;
		real_t *dB = nullptr;

		virtual const real_t *EvaluateDeltaBOverB(const real_t) override;

	public:
		HeatTransportRRAdaptiveMHDLike(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			const real_t, const real_t, const real_t
		);
		virtual ~HeatTransportRRAdaptiveMHDLike();
	};
}

#endif/*_DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_ADAPTIVE_MHD_LIKE_HPP*/
