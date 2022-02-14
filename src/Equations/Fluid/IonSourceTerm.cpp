/**
 * Implementation of an equation term allowing for prescribed ion densities.
 */

#include "DREAM/Equations/Fluid/IonSourceTerm.hpp"
#include "DREAM/MultiInterpolator1D.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
IonSourceTerm::IonSourceTerm(
    FVM::Grid *grid, IonHandler *ihdl, const len_t nIons,
    const len_t *ionIndices, MultiInterpolator1D *data
) : IonPrescribedParameter(grid, ihdl, nIons, ionIndices, data) {

    SetName("IonSourceTerm");
}


/**
 * Destructor.
 */
IonSourceTerm::~IonSourceTerm() {}


/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 * x:       Current value of the unknown quantity.
 *
 * (This term represents a constant, and since the derivative
 * with respect to anything of a constant is zero, we don't need
 * to do anything).
 */
bool IonSourceTerm::SetJacobianBlock(
    const len_t, const len_t, FVM::Matrix*, const real_t*
) {
	// Derivative of constant = 0
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
void IonSourceTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    const len_t Nr = this->grid->GetNr();

    for (len_t i = 0; i < nIons; i++) {
        for (len_t Z0 = 0; Z0 <= Z[i]; Z0++) {
            const len_t idx = this->ions->GetIndex(ionIndices[i], Z0);
            real_t *n = currentData[i] + Z0*Nr;

            for (len_t ir = 0; ir < Nr; ir++)
                rhs[idx*Nr+ir] += n[ir];
        }
    }
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * ni:  Ion densities in previous iteration.
 */
void IonSourceTerm::SetVectorElements(real_t *vec, const real_t*) {
    const len_t Nr = this->grid->GetNr();

    for (len_t i = 0; i < nIons; i++) {
        for (len_t Z0 = 0; Z0 <= Z[i]; Z0++) {
            const len_t idx = this->ions->GetIndex(ionIndices[i], Z0);
            real_t *n = currentData[i] + Z0*Nr;

            for (len_t ir = 0; ir < Nr; ir++)
                vec[idx*Nr+ir] += n[ir];
        }
    }
}

