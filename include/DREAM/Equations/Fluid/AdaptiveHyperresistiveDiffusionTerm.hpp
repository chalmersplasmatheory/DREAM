#ifndef _DREAM_EQUATIONS_FLUID_ADAPTIVE_HYPERRESISTIVE_DIFFUSION_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_ADAPTIVE_HYPERRESISTIVE_DIFFUSION_TERM_HPP

#include "DREAM/Equations/AdaptiveMHDLikeTransportTerm.hpp"
#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"

namespace DREAM {
	class AdaptiveHyperresistiveDiffusionTerm
		: public AdaptiveMHDLikeTransportTerm,
		  public HyperresistiveDiffusionTerm {
	
	protected:
		real_t Lambda0;
		real_t *Lambda;

		virtual const real_t *EvaluateLambda(const real_t t) override;

	public:
		AdaptiveHyperresistiveDiffusionTerm(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			const real_t, bool, const real_t, const real_t
		);
		virtual ~AdaptiveHyperresistiveDiffusionTerm();

		virtual const real_t *GetLambda() const override { return this->Lambda; }
	};
}

#endif/*_DREAM_EQUATIONS_FLUID_ADAPTIVE_HYPERRESISTIVE_DIFFUSION_TERM_HPP*/
