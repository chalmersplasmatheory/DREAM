#ifndef _DREAM_EQUATIONS_ADAPTIVE_MHD_LIKE_TRANSPORT_TERM_HPP
#define _DREAM_EQUATIONS_ADAPTIVE_MHD_LIKE_TRANSPORT_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"


namespace DREAM {
	class AdaptiveMHDLikeTransportTerm {
	protected:
		FVM::Grid *grid;
		FVM::UnknownQuantityHandler *uqh;

		// Maximum current density gradient
		real_t grad_j_tot_max;
		bool gradient_normalized = false;

		bool transport_enabled = false;
		real_t transport_enabled_t = 0;
		real_t min_duration = 0;

		len_t id_j_tot;

		real_t volume = 0;
		real_t javg = 0;

		bool IsCurrentGradientExceeded(const len_t);
	public:
		AdaptiveMHDLikeTransportTerm(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			const real_t, bool, const real_t
		);

		bool CheckTransportEnabled(const real_t);

		bool IsCurrentGradientExceeded();
	};
}

#endif/*_DREAM_EQUATIONS_ADAPTIVE_MHD_LIKE_TRANSPORT_TERM_HPP*/
