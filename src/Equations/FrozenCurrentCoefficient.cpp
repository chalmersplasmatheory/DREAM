/**
 * Implementation of the equation for the frozen current coefficient.
 */

#include "DREAM/Equations/FrozenCurrentCoefficient.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Matrix.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
FrozenCurrentCoefficient::FrozenCurrentCoefficient(
	FVM::Grid *grid, FVM::Grid *fluidGrid, FVM::Interpolator1D *I_p_presc,
	FVM::UnknownQuantityHandler *unknowns, const real_t D_I_max
) : EquationTerm(grid), fluidGrid(fluidGrid), I_p_presc(I_p_presc),
	id_I_p(unknowns->GetUnknownID(OptionConstants::UQTY_I_P)),
	id_j_tot(unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT)),
	D_I_max(D_I_max) { }


/**
 * Destructor.
 */
FrozenCurrentCoefficient::~FrozenCurrentCoefficient() {
	if (this->I_p_presc != nullptr)
		delete this->I_p_presc;
}

/**
 * Rebuild this equation term.
 */
void FrozenCurrentCoefficient::Rebuild(
	const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns
) {
	real_t Ipresc = this->I_p_presc->Eval(t)[0];
	real_t Ip = unknowns->GetUnknownData(this->id_I_p)[0];
	real_t DIk = this->D_I;

	real_t Irat;
	if (Ipresc != 0)
		Irat = Ip/Ipresc;
	else
		Irat = 2;

	if (this->D_I == 0 && (Irat <= 1)) {
		this->D_I = 0;
	} else if (this->D_I == 0) {
		// Estimate D_I based on formula
		const len_t nr = this->fluidGrid->GetNr();
		const real_t *dr = this->fluidGrid->GetRadialGrid()->GetDr();

		real_t Ip_old = unknowns->GetUnknownDataPrevious(this->id_I_p)[0];
		real_t dIpdt = (Ip - Ip_old) / dt;
		real_t VpVol = this->fluidGrid->GetVpVol(nr-1);
		real_t ja = unknowns->GetUnknownData(this->id_j_tot)[nr-1];

		// assume j(r>a) = 0
		real_t djdr = ja / dr[nr-1];

		this->D_I = 2*M_PI / VpVol * dIpdt / djdr;
	} else if (this->prevTime != t) {
		this->D_I *= Irat;
	} else {
		this->D_I -= (DIk - this->D_I_prev) * Ipresc * (Irat - 1) / (Ip - Ip_prev);
	}

	// Bound lower value to zero
	if (this->D_I < 0)
		this->D_I = 0;
	// Bound upper value
	if (this->D_I > this->D_I_max)
		this->D_I = this->D_I_max;
	
	this->D_I_prev = DIk;
	this->Ip_prev = Ip;
	this->prevTime = t;
}

/**
 * Set the jacobian elements for this equation term.
 */
bool FrozenCurrentCoefficient::SetJacobianBlock(
	const len_t, const len_t, FVM::Matrix*, const real_t*
) {
	throw FVM::EquationTermException(
		"FrozenCurrentCoefficient: This term should be added to the external "
		"iterator."
	);
}

/**
 * Set the matrix elements for this equation term.
 */
void FrozenCurrentCoefficient::SetMatrixElements(
	FVM::Matrix*, real_t*
) {
	throw FVM::EquationTermException(
		"FrozenCurrentCoefficient: This term should be added to the external "
		"iterator."
	);
}

/**
 * Set the non-linear vector elements for this equation term.
 */
void FrozenCurrentCoefficient::SetVectorElements(
	real_t *F, const real_t*
) {
	F[0] = this->D_I;
}

