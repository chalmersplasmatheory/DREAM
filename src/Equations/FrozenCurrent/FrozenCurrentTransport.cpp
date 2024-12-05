/**
 * Implementation of the transport operator used for the "frozen current" mode,
 * where the radial transport rate is gradually adjusted so that the simulated
 * total plasma current meaures a prescribed value.
 */

#include "DREAM/Equations/FrozenCurrentTransport.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
FrozenCurrentTransport::FrozenCurrentTransport(
	FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns,
	enum TransportMode ts, bool allocCoefficients
) : DiffusionTerm(g, allocCoefficients), transportMode(ts),
	id_D_I(unknowns->GetUnknownID(OptionConstants::UQTY_D_I)) {

	AddUnknownForJacobian(unknowns, id_D_I);
}


/**
 * Destructor.
 */
FrozenCurrentTransport::~FrozenCurrentTransport() {}

/**
 * Method which is called when the simulation grid is rebuilt.
 */
bool FrozenCurrentTransport::GridRebuilt() {
	this->DiffusionTerm::GridRebuilt();

	return true;
}

/**
 * Method for rebuilding the diffusion coefficient of this operator.
 */
void FrozenCurrentTransport::Rebuild(
	const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler *unknowns
) {
	const real_t D_I = unknowns->GetUnknownData(this->id_D_I)[0];

	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	const real_t
		*p   = mg->GetP1(),
		*xi0 = mg->GetP2();
	const len_t
		np  = mg->GetNp1(),
		nxi = mg->GetNp2();
	
	const bool kinetic = (np > 1 || nxi > 1);

	for (len_t ir = 0; ir < nr+1; ir++) {
		for (len_t i = 0; i < np; i++) {
			if (this->transportMode == TRANSPORT_MODE_CONSTANT) {
				for (len_t j = 0; j < nxi; j++)
					Drr(ir,i,j) += D_I;
			} else if (this->transportMode == TRANSPORT_MODE_BETAPAR) {
				if (kinetic) {
					real_t beta = p[i] / sqrt(1+p[i]*p[i]);
					for (len_t j = 0; j < nxi; j++)
						Drr(ir,i,j) += D_I * beta * xi0[j];
				} else {
					Drr(ir,0,0) += D_I;
				}
			}
		}
	}
}

/**
 * Set the partial derivative of the diffusion coefficient.
 *
 * derivId:    ID of the unknown quantity to differentiate w.r.t.
 * nMultiples: Number of multiples in unknown quantity (should always be 1
 *             for this term).
 */
void FrozenCurrentTransport::SetPartialDiffusionTerm(
	len_t derivId, len_t
) {
	if (derivId != id_D_I)
		return;

	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	const real_t
		*p   = mg->GetP1(),
		*xi0 = mg->GetP2();
	const len_t
		np  = mg->GetNp1(),
		nxi = mg->GetNp2();

	const bool kinetic = (np > 1 || nxi > 1);

	for (len_t ir = 0; ir < nr+1; ir++) {
		for (len_t i = 0; i < np; i++) {
			if (this->transportMode == TRANSPORT_MODE_CONSTANT) {
				for (len_t j = 0; j < nxi; j++)
					dDrr(ir,i,j,0) = 1;
			} else if (this->transportMode == TRANSPORT_MODE_BETAPAR) {
				if (kinetic) {
					real_t beta = p[i] / sqrt(1+p[i]*p[i]);
					for (len_t j = 0; j < nxi; j++)
						dDrr(ir,i,j,0) = beta * xi0[j];
				} else {
					dDrr(ir,0,0,0) = 1;
				}
			}
		}
	}
}

