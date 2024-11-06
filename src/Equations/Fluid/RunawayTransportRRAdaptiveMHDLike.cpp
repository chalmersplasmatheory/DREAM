/**
 */

#include "DREAM/Equations/Fluid/RunawayTransportRRAdaptiveMHDLike.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
RunawayTransportRRAdaptiveMHDLike::RunawayTransportRRAdaptiveMHDLike(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *uqh,
	const real_t grad_j_tot_max, const real_t min_duration,
	const real_t dBB0
) : AdaptiveMHDLikeTransportTerm(grid, uqh, grad_j_tot_max, min_duration),
	RunawayTransportRechesterRosenbluth(grid, nullptr), dBOverB(dBB0) {
		
	this->dB = new real_t[grid->GetNr()];
}


/**
 * Destructor.
 */
RunawayTransportRRAdaptiveMHDLike::~RunawayTransportRRAdaptiveMHDLike() {
	if (this->dB != nullptr)
		delete [] this->dB;
}


/**
 * Evaluate dB/B.
 */
const real_t *RunawayTransportRRAdaptiveMHDLike::EvaluateDeltaBOverB(const real_t t) {
	real_t v = 0;
	if (this->CheckTransportEnabled(t))
		v = this->dBOverB;
	
	const len_t nr = this->AdaptiveMHDLikeTransportTerm::grid->GetNr();
	for (len_t ir = 0; ir < nr; ir++)
		this->dB[ir] = v;
	
	return this->dB;
}


