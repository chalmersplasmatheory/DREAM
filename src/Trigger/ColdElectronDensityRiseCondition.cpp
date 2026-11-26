
#include "DREAM/Trigger/ColdElectronDensityRiseCondition.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;


/**
 * Constructor.
 *
 * sensitivity: The condition is triggered when n_cold/n_hot > sensitivity.
 */
ColdElectronDensityRiseCondition::ColdElectronDensityRiseCondition(
	FVM::Grid *g, FVM::UnknownQuantityHandler *u,
	const len_t nMultiples, const real_t sensitivity,
	bool saveTriggerState
) : EquationTriggerCondition(g, u, nMultiples, saveTriggerState), sensitivity(sensitivity) {

	this->id_n_cold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->id_n_hot = u->GetUnknownID(OptionConstants::UQTY_N_HOT);
}

/**
 * Destructor.
 */
ColdElectronDensityRiseCondition::~ColdElectronDensityRiseCondition() {
}


/**
 * Check whether the trigger condition is enabled.
 */
void ColdElectronDensityRiseCondition::CheckCondition(
	const real_t, FVM::UnknownQuantityHandler *unknowns
) {
	const real_t *n_cold = unknowns->GetUnknownDataPrevious(this->id_n_cold);
	const real_t *n_hot = unknowns->GetUnknownDataPrevious(this->id_n_hot);

	const len_t nr = this->grid->GetNr();
	const len_t Np = this->grid->GetNp1(0) * this->grid->GetNp2(0);
	for (len_t in = 0, i = 0; in < this->nMultiples; in++) {
		for (len_t ir = 0; ir < nr; ir++) {
			bool v;
			if (n_hot[ir] < 1)
				v = true;
			else {
				real_t r = n_cold[ir] / n_hot[ir];
				v = (r > sensitivity);
			}

			for (len_t ip = 0; ip < Np; ip++, i++)
				triggered[i] = v;
		}
	}
}

