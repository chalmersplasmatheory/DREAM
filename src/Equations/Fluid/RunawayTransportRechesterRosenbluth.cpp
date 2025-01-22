/**
 * Implementation of a Rechester-Rosenbluth diffusion operator for the
 * runaway electron density, n_re.
 */

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/RunawayTransportRechesterRosenbluth.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
RunawayTransportRechesterRosenbluth::RunawayTransportRechesterRosenbluth(
	FVM::Grid *grid, FVM::Interpolator1D *dBB
) : FVM::DiffusionTerm(grid), dBB(dBB) { }


/**
 * Destructor.
 */
RunawayTransportRechesterRosenbluth::~RunawayTransportRechesterRosenbluth() {
	if (this->dBB != nullptr)
		delete this->dBB;
}


/**
 * Rebuild the diffusion coefficient for this term.
 */
void RunawayTransportRechesterRosenbluth::Rebuild(
	const real_t t, const real_t, FVM::UnknownQuantityHandler*
) {
	const real_t *dB_B = this->EvaluateDeltaBOverB(t);
	FVM::RadialGrid *rg = this->grid->GetRadialGrid();
	const len_t nr = this->grid->GetNr();

	real_t qR0;
	const real_t R0 = rg->GetR0();
	if (isinf(R0))
		qR0 = 1;
	else
		qR0 = R0;

	for (len_t ir = 0; ir < nr; ir++)
		Drr(ir, 0, 0) += M_PI * qR0 * Constants::c * dB_B[ir];
}


