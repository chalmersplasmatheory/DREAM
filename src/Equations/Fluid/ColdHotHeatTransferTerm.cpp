/**
 * An equation term representing the transfer of heat from the cold to the hot
 * region in fully kinetic mode when the cold/hot boundary definition changes.
 */

#include "DREAM/Equations/Fluid/ColdHotHeatTransferTerm.hpp"
#include "DREAM/Constants.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
ColdHotHeatTransferTerm::ColdHotHeatTransferTerm(
	FVM::Grid *fluidGrid, FVM::Grid *kineticGrid,
	len_t id_W, len_t id_f,
	FVM::UnknownQuantityHandler *u, real_t pThreshold,
	pThresholdMode pMode, real_t scaleFactor
) : MomentQuantity(fluidGrid, kineticGrid, id_W, id_f, u, pThreshold, pMode),
	scaleFactor(scaleFactor) {

	const len_t nr = kineticGrid->GetNr();
	const len_t np = kineticGrid->GetMomentumGrid(0)->GetNp1();
	this->prevThresholdEnvelope = new real_t[nr*np];
	this->prevThresholdValue = new real_t[nr];

	StoreThresholdEnvelope(false);
}

/**
 * Destructor.
 */
ColdHotHeatTransferTerm::~ColdHotHeatTransferTerm() {
	delete [] this->prevThresholdValue;
	delete [] this->prevThresholdEnvelope;
}


/**
 * Method that is called whenever the grid is rebuilt.
 */
bool ColdHotHeatTransferTerm::GridRebuilt() {
	if (this->MomentQuantity::GridRebuilt()) {
		delete [] this->prevThresholdValue;
		delete [] this->prevThresholdEnvelope;

		// XXX Here we assume that momentum grids are the same at all radii
		const len_t nr = fGrid->GetNr();
		const len_t np = fGrid->GetMomentumGrid(0)->GetNp1();
		this->prevThresholdEnvelope = new real_t[nr*np];
		this->prevThresholdValue = new real_t[nr];
		StoreThresholdEnvelope(false);

		return true;
	} else
		return false;
}

void ColdHotHeatTransferTerm::Rebuild(
	const real_t, const real_t, FVM::UnknownQuantityHandler*
) {
	// XXX: Here we assume that momentum grids are the same at all radii
	const len_t nr = fGrid->GetNr();
	const len_t np1 = fGrid->GetMomentumGrid(0)->GetNp1();
	const len_t np2 = fGrid->GetMomentumGrid(0)->GetNp2();

	const real_t mc2 = Constants::me*Constants::c*Constants::c;

	len_t offset = 0;
	for (len_t ir = 0; ir < nr; ir++) {
		FVM::MomentumGrid *mg = fGrid->GetMomentumGrid(ir);
		const real_t *env = this->prevThresholdEnvelope+(ir*np1);

		for (len_t i = 0; i < np1; i++) {
			for (len_t j = 0; j < np2; j++) {
				len_t ind = j*np1 + i;
				real_t g = mg->GetGamma(i, j);
				this->integrand[offset+ind] =
					this->scaleFactor * (1-env[i]) * mc2 * (g-1);
			}
		}

		offset += np1*np2;
	}

	// Evaluate threshold envelope in this time step
	StoreThresholdEnvelope();
}

/**
 * Evaluate and store the threshold envelope value at all momentum
 * grid points (p only, not xi).
 */
void ColdHotHeatTransferTerm::StoreThresholdEnvelope(bool updateIntegrationDirection) {
	const len_t nr = fGrid->GetNr();
	const len_t np = fGrid->GetMomentumGrid(0)->GetNp1();
	for (len_t ir = 0; ir < nr; ir++) {
		for (len_t i = 0; i < np; i++) {
			this->prevThresholdEnvelope[ir*np+i] = this->ThresholdEnvelope(ir, i);
		}

		real_t p0 = this->ThresholdValue(ir);

		if (updateIntegrationDirection) {
			if (this->prevThresholdValue[ir] > p0)
				this->pIntMode[ir] = P_INT_MODE_THR2MAX;
			else
				this->pIntMode[ir] = P_INT_MODE_ZER2THR;
		}

		this->prevThresholdValue[ir] = p0;
	}
}

