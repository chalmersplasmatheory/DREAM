/**
 * Adaptive MHD-like heat transport, using a Rechester-Rosenbluth
 * diffusion coefficient.
 */

#include "DREAM/Equations/Fluid/HeatTransportRRAdaptiveMHDLike.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
HeatTransportRRAdaptiveMHDLike::HeatTransportRRAdaptiveMHDLike(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *uqh,
	const real_t grad_j_tot_max, bool gradient_normalized,
	const real_t dBOverB, const real_t suppression_level, bool localized
) : AdaptiveMHDLikeTransportTerm(grid, uqh, grad_j_tot_max, gradient_normalized, suppression_level, localized),
	HeatTransportRechesterRosenbluth(grid, nullptr, uqh),
	dBOverB(dBOverB) {

	this->dB = new real_t[grid->GetNr()];
}


/**
 * Destructor.
 */
HeatTransportRRAdaptiveMHDLike::~HeatTransportRRAdaptiveMHDLike() {
	if (this->dB != nullptr)
		delete [] this->dB;
}


/**
 * Evaluate dB/B.
 */
const real_t *HeatTransportRRAdaptiveMHDLike::EvaluateDeltaBOverB(const real_t t) {
	real_t v = 0;
	if (this->CheckTransportEnabled(t))
		v = this->dBOverB;
	
	const len_t nr = this->AdaptiveMHDLikeTransportTerm::grid->GetNr();
	for (len_t ir = 0; ir < nr; ir++)
		this->dB[ir] = v * this->mask[ir];
	
	return this->dB;
}


