/**
 * Implementation of the equation for the frozen current coefficient.
 */

#include "DREAM/Equations/FrozenCurrentCoefficient_D.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
FrozenCurrentCoefficient_D::FrozenCurrentCoefficient_D(
	FVM::Grid *grid, FVM::Interpolator1D *I_p_presc,
	FVM::UnknownQuantityHandler *unknowns
) : EquationTerm(grid), I_p_presc(I_p_presc),
	id_D_I(unknowns->GetUnknownID(OptionConstants::UQTY_D_I)),
	id_I_p(unknowns->GetUnknownID(OptionConstants::UQTY_I_P)) { }


/**
 * Destructor.
 */
FrozenCurrentCoefficient_D::~FrozenCurrentCoefficient_D() {}

/**
 * Rebuild this equation term.
 */
void FrozenCurrentCoefficient_D::Rebuild(
	const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
	// If this is the first iteration of the time step, check whether to go
	// to prescribed current mode.
	if (!this->prescribeCurrent && t > this->prevTime) {
		real_t Ipresc = this->I_p_presc->Eval(t)[0];
		real_t Ip = unknowns->GetUnknownData(id_I_p)[0];
		if (Ip >= Ipresc)
			this->prescribeCurrent = true;
		
		this->prevTime = t;
	}
}

/**
 * Set the jacobian elements for this equation term.
 */
bool FrozenCurrentCoefficient_D::SetJacobianBlock(
	const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*
) {
	// D_I = 0
	if (!this->prescribeCurrent && derivId == id_D_I) {
		jac->SetElement(0, 0, 1);
		return true;
	}

	return false;
}

/**
 * Set the matrix elements for this equation term.
 */
void FrozenCurrentCoefficient_D::SetMatrixElements(
	FVM::Matrix *mat, real_t *rhs
) {
	if (!this->prescribeCurrent) {
		mat->SetElement(0, 0, 1.0);
		rhs[0] = 0;
	}
}

/**
 * Set the non-linear vector elements for this equation term.
 */
void FrozenCurrentCoefficient_D::SetVectorElements(
	real_t *F, const real_t *D_I
) {
	if (!this->prescribeCurrent) {
		// D_I = 0
		F[0] = *D_I;
	}
}

