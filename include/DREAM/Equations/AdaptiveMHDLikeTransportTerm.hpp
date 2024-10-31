#ifndef _DREAM_EQUATIONS_ADAPTIVE_MHD_LIKE_TRANSPORT_TERM_HPP
#define _DREAM_EQUATIONS_ADAPTIVE_MHD_LIKE_TRANSPORT_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"


namespace DREAM {
	class AdaptiveMHDLikeTransportTerm : public FVM::EquationTerm {
	protected:
		FVM::UnknownQuantityHandler *uqh;

		// Maximum current density gradient
		real_t grad_j_tot_max;

		bool transport_enabled = false;
		real_t transport_enabled_t = 0;
		real_t min_duration = 0;

		len_t id_j_tot;

	public:
		AdaptiveMHDLikeTransportTerm(
			FVM::Grid*, FVM::UnknownQuantityHandler*, const real_t,
			const real_t
		);

		bool CheckTransportEnabled(const real_t);

		bool IsCurrentGradientExceeded();
		bool IsCurrentGradientExceeded(const len_t);
	};
}

#endif/*_DREAM_EQUATIONS_ADAPTIVE_MHD_LIKE_TRANSPORT_TERM_HPP*/
