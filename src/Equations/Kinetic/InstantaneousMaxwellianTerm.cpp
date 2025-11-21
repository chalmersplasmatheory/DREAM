/**
 * Implementation of a kinetic term which gives a Maxwellian evaluated at the
 * instantaneous density and temperature in every time step.
 */

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/InstantaneousMaxwellianTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
InstantaneousMaxwellianTerm::InstantaneousMaxwellianTerm(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *uqh,
	enum MaxwellianPopulation pop
) : EquationTerm(grid) {

	SetName("InstantaneousMaxwellianTerm");

	if (pop == MAXWELLIAN_POPULATION_HOT) {
		this->id_T = uqh->GetUnknownID(OptionConstants::UQTY_T_HOT);
		this->id_n = uqh->GetUnknownID(OptionConstants::UQTY_N_HOT);
	} else {
		this->id_T = uqh->GetUnknownID(OptionConstants::UQTY_T_COLD);
		this->id_n = uqh->GetUnknownID(OptionConstants::UQTY_N_COLD);
	}

    this->F = new real_t[grid->GetNCells()];
	this->dFdn = new real_t[grid->GetNCells()];
	this->dFdT = new real_t[grid->GetNCells()];
}

/**
 * Destructor.
 */
InstantaneousMaxwellianTerm::~InstantaneousMaxwellianTerm() {
	delete [] this->F;
	delete [] this->dFdn;
	delete [] this->dFdT;
}


/**
 * Build the value for this quantity in the next time step.
 */
void InstantaneousMaxwellianTerm::Rebuild(
	const real_t, const real_t, FVM::UnknownQuantityHandler *uqh
) {
	const real_t *T = uqh->GetUnknownData(this->id_T);
	const real_t *n = uqh->GetUnknownData(this->id_n);

	const len_t nr = this->grid->GetNr();
	const len_t np = this->grid->GetNp1(0);
	const len_t nxi = this->grid->GetNp2(0);
	const real_t *p = this->grid->GetMomentumGrid(0)->GetP1();

	for (len_t ir = 0; ir < nr; ir++) {
		for (len_t ip = 0; ip < np; ip++) {
			const len_t idx0 = ir*np*nxi + ip;
			this->F[idx0] = Constants::RelativisticMaxwellian(
				p[ip], n[ir], T[ir],
				this->dFdn+idx0, this->dFdT+idx0
			);

			// Copy to all other xi
			for (len_t ix = 1; ix < nxi; ix++) {
				this->F[idx0 + ix*np] = this->F[idx0];
				this->dFdn[idx0 + ix*np] = this->dFdn[idx0];
				this->dFdT[idx0 + ix*np] = this->dFdT[idx0];
			}
		}
	}
}


/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * jac:     Jacobian matrix block to populate.
 * x:       Value of the unknown quantity.
 */
bool InstantaneousMaxwellianTerm::SetJacobianBlock(
	const len_t, const len_t derivId,
	FVM::Matrix *jac, const real_t*
) {
	const real_t *dd = nullptr;
	if (derivId == this->id_n)
		dd = this->dFdn;
	else if (derivId == this->id_T)
		dd = this->dFdT;
	
	if (dd != nullptr) {
		const len_t nr = this->grid->GetNr();
		const len_t np = this->grid->GetNp1(0);
		const len_t nxi = this->grid->GetNp2(0);

		for (len_t ir = 0, i = 0; ir < nr; ir++)
			for (len_t j = 0; j < nxi*np; j++, i++)
				jac->SetElement(i, ir, dd[i]);

		return true;
	} else
		return false;
}


/**
 * Set the elements in the matrix and on the RHS corresponding
 * to this quantity.
 *
 * mat: Matrix to set elements in (1 is added to the diagonal)
 * rhs: Right-hand-side. Values will be set to the current value of
 *      this parameter.
 */
void InstantaneousMaxwellianTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        rhs[i] = -this->F[i];
}


/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution.
 */
void InstantaneousMaxwellianTerm::SetVectorElements(real_t *vec, const real_t*) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] = this->F[i];
}

