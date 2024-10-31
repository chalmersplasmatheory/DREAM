/**
 * Implementation of the hyperresistive diffusion term with an adaptive
 * diffusion coefficient.
 */

#include "DREAM/Equations/Fluid/AdaptiveHyperresistiveDiffusionTerm.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
AdaptiveHyperresistiveDiffusionTerm::AdaptiveHyperresistiveDiffusionTerm(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *uqh,
	const real_t grad_j_tot_max, const real_t Lambda0,
	const real_t min_duration
) : AdaptiveMHDLikeTransportTerm(grid, uqh, grad_j_tot_max, min_duration),
	HyperresistiveDiffusionTerm(grid, nullptr), Lambda0(Lambda0) {
	
	this->Lambda = new real_t[grid->GetNr()];
}


/**
 * Destructor.
 */
AdaptiveHyperresistiveDiffusionTerm::~AdaptiveHyperresistiveDiffusionTerm() {
	delete [] this->Lambda;
}


/**
 * Return the diffusion coefficient.
 *
 * t: Time to evaluate the coefficient at.
 */
const real_t *AdaptiveHyperresistiveDiffusionTerm::EvaluateLambda(const real_t t) {
	real_t v = 0;
	if (this->CheckTransportEnabled(t))
		v = this->Lambda0;
	
	const len_t nr = this->AdaptiveMHDLikeTransportTerm::grid->GetNr();
	for (len_t ir = 0; ir < nr; ir++)
		this->Lambda[ir] = v;
	
	return this->Lambda;
}


