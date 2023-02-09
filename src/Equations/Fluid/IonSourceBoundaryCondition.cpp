/**
 * Implementation of an ion source boundary condition.
 */

#include "DREAM/Equations/Fluid/IonSourceBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryCondition.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
IonSourceBoundaryCondition::IonSourceBoundaryCondition(
	FVM::Grid *grid, IonHandler *ions, MultiInterpolator1D *dndt,
	const len_t iIon
) : IonBoundaryCondition<FVM::BC::BoundaryCondition>(grid, ions, iIon), dndt(dndt) {
}

/**
 * Destructor.
 */
IonSourceBoundaryCondition::~IonSourceBoundaryCondition() {
	delete this->dndt;
}

/**
 * Rebuild this boundary conditon.
 */
bool IonSourceBoundaryCondition::Rebuild(
	const real_t t, FVM::UnknownQuantityHandler*
) {
	this->time = t;

	return true;
}


/**
 * Add elements to linearly implicit equation system.
 */
void IonSourceBoundaryCondition::AddCSToMatrixElements(
	FVM::Matrix*, real_t *rhs, const len_t iIon,
	const len_t Z0, const len_t rOffset
) {
	const len_t idx = this->ions->GetIndex(iIon, Z0);
	const len_t Nr = this->grid->GetNr();
	const real_t *dn = dndt->Eval(idx, this->time);
	const real_t VpVol = this->grid->GetVpVol()[Nr-1];
	const real_t dr = this->grid->GetRadialGrid()->GetDr(Nr-1);

	rhs[rOffset + (Nr-1)] += dn[0] / (VpVol*dr);
}

/**
 * Add elements to non-linear function vector.
 */
void IonSourceBoundaryCondition::AddCSToVectorElements(
	real_t *vec, const real_t*, const len_t iIon,
	const len_t Z0, const len_t rOffset
) {
	const len_t idx = this->ions->GetIndex(iIon, Z0);
	const real_t *dn = dndt->Eval(idx, this->time);
	const len_t Nr = this->grid->GetNr();
	const real_t VpVol = this->grid->GetVpVol()[Nr-1];
	const real_t dr = this->grid->GetRadialGrid()->GetDr(Nr-1);

	vec[rOffset + (Nr-1)] += dn[0] / (VpVol*dr);
}

