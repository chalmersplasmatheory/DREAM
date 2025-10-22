
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
	const real_t sensitivity
) : EquationTriggerCondition(g, u), sensitivity(sensitivity) {

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
	const real_t *n_cold = unknowns->GetUnknownData(this->id_n_cold);
	const real_t *n_hot = unknowns->GetUnknownData(this->id_n_hot);

	const len_t N = this->grid->GetNCells();
	for (len_t i = 0; i < N; i++) {
		if (n_hot[i] < 1)
			triggered[i] = true;
		else {
			real_t r = n_cold[i] / n_hot[i];
			triggered[i] = (r > sensitivity);
		}
	}
}

