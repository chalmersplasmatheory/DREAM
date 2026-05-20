#ifndef _DREAM_EQUATIONS_FLUID_ADAPTIVE_HYPERRESISTIVE_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_ADAPTIVE_HYPERRESISTIVE_DIFFUSION_TERM_HPP

#include "DREAM/Equations/AdaptiveMHDLikeTransportTerm.hpp"
#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"

namespace DREAM {
	class AdaptiveHyperresistiveDiffusionTerm
		: public AdaptiveMHDLikeTransportTerm,
		  public HyperresistiveDiffusionTerm {
	
	protected:
		IonHandler *ions;
		real_t dBB0;
		real_t *Lambda, *Lambda0, *dLambda;

		len_t id_n_i;

		virtual const real_t *EvaluateLambda(const real_t t) override;

	public:
		AdaptiveHyperresistiveDiffusionTerm(
			FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*,
			const real_t, bool, const real_t, const real_t, bool
		);
		virtual ~AdaptiveHyperresistiveDiffusionTerm();

		virtual const real_t *GetLambda() const override { return this->Lambda; }
		virtual void SetPartialDiffusionTerm(len_t, len_t) override;
	};
}

#endif/*_DREAM_EQUATIONS_FLUID_ADAPTIVE_HYPERRESISTIVE_DIFFUSION_TERM_HPP*/
