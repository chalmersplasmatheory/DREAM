#ifndef _DREAM_EQUATIONS_RUNAWAY_TRANSPORT_RR_ADAPTIVE_MHD_LIKE_HPP
#define _DREAM_EQUATIONS_RUNAWAY_TRANSPORT_RR_ADAPTIVE_MHD_LIKE_HPP

#include "DREAM/Equations/AdaptiveMHDLikeTransportTerm.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/Equations/Fluid/RunawayTransportRechesterRosenbluth.hpp"

namespace DREAM {
	class RunawayTransportRRAdaptiveMHDLike :
		public AdaptiveMHDLikeTransportTerm,
		public RunawayTransportRechesterRosenbluth {
	protected:
		real_t dBOverB;
		real_t *dB = nullptr;

		virtual const real_t *EvaluateDeltaBOverB(const real_t) override;

	public:
		RunawayTransportRRAdaptiveMHDLike(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			const real_t, const real_t, const real_t
		);
		virtual ~RunawayTransportRRAdaptiveMHDLike();
	};
}

#endif/*_DREAM_EQUATIONS_RUNAWAY_TRANSPORT_RR_ADAPTIVE_MHD_LIKE_HPP*/
