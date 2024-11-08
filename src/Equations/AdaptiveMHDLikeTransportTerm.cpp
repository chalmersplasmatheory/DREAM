/**
 * The 'AdaptiveMHDLikeTransportTerm' is an abstract class which can be used
 * for the implementation of transport processes which are induced by MHD.
 */

#include <cmath>
#include "DREAM/Equations/AdaptiveMHDLikeTransportTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
AdaptiveMHDLikeTransportTerm::AdaptiveMHDLikeTransportTerm(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *uqh,
	const real_t grad_j_tot_max, bool gradient_normalized,
	const real_t min_duration, bool localized
) : grid(grid), uqh(uqh),
	grad_j_tot_max(grad_j_tot_max),
	gradient_normalized(gradient_normalized),
	min_duration(min_duration), localized(localized),
	id_j_tot(uqh->GetUnknownID(OptionConstants::UQTY_J_TOT)) {
	
	this->mask = new real_t[grid->GetNr()];
}


/*
 * Destructor.
 */
AdaptiveMHDLikeTransportTerm::~AdaptiveMHDLikeTransportTerm() {
	delete [] this->mask;
}


/**
 * Check if the transport should be enabled. This routine returns
 * true from when the current density gradient exceeds the threshold
 * value, until the gradient has been sufficiently reduced, or for
 * at least 'min_duration' seconds.
 */
bool AdaptiveMHDLikeTransportTerm::CheckTransportEnabled(const real_t t) {
	if (this->gradient_normalized) {
		const real_t *j_tot = this->uqh->GetUnknownDataPrevious(this->id_j_tot);
		const real_t *dr = grid->GetRadialGrid()->GetDr();
		const real_t a = grid->GetRadialGrid()->GetMinorRadius();
		this->javg = 0;
		for (len_t ir = 0; ir < grid->GetNr(); ir++)
			this->javg += j_tot[ir] * dr[ir];

		this->javg /= a;
	}

	if (this->transport_enabled) {
		// Disable transport?
		// (we disable transport if 'min_duration' has been exceeded,
		// AND the current density gradient remains high)
		real_t dt = t - this->transport_enabled_t;
		if (dt > this->min_duration)
			this->transport_enabled = IsCurrentGradientExceeded();
	} else {
		// Enable transport
		this->transport_enabled = IsCurrentGradientExceeded();
		this->transport_enabled_t = t;
	}

	return this->transport_enabled;
}


/**
 * Checks whether the current density gradient is
 * exceeded anywhere on the grid.
 */
bool AdaptiveMHDLikeTransportTerm::IsCurrentGradientExceeded() {
	const len_t nr = this->grid->GetNr();
	bool exceeded = false;

	for (len_t ir = 0; ir < nr && !exceeded; ir++) {
		bool e = this->IsCurrentGradientExceeded(ir);

		if (this->localized) {
			if (e) {
				this->mask[ir] = 1;
				if (ir > 0)
					this->mask[ir-1] = 1;
				if (ir < nr-1)
					this->mask[ir+1] = 1;

				exceeded = true;
			} else
				this->mask[ir] = 0;
		} else
			// Apply uniformly everywhere
			this->mask[ir] = 1;
	}
	
	return exceeded;
}


/**
 * Checks whether the current density gradient is
 * exceeded in the specified radial grid point.
 *
 * ir: Radial grid point to check current density gradient in.
 */
bool AdaptiveMHDLikeTransportTerm::IsCurrentGradientExceeded(const len_t ir) {
	const real_t *j_tot = this->uqh->GetUnknownDataPrevious(this->id_j_tot);
	const real_t *dr_f = this->grid->GetRadialGrid()->GetDr_f();

	// Calculate backward derivatives, except at ir = 0
	real_t gradj;
	if (ir == 0)
		gradj = (j_tot[1]-j_tot[0]) / dr_f[0];
	else
		gradj = (j_tot[ir]-j_tot[ir-1]) / dr_f[ir-1];
	
	if (this->gradient_normalized) {
		const real_t a = this->grid->GetRadialGrid()->GetMinorRadius();
		gradj *= a / this->javg;
	}
	
	return (std::abs(gradj) >= this->grad_j_tot_max);
}

