/**
 * Implementation of a general ion equation term, which applies to
 * all ions and ion charge states in the equation system.
 */

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
IonEquationTerm::IonEquationTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon)
    : FVM::EquationTerm(g), ions(ihdl), iIon(iIon) {

    this->Zion = this->ions->GetZ(iIon);
}

/**
 * Destructor.
 */
IonEquationTerm::~IonEquationTerm() {}


/**
 * Set the specified block in the given Jacobian matrix.
 *
 * derivId: ID of the unknown quantity with respect to
 *          which differentiation should be done.
 * uqtyId:  ID of the unknown quantity to differentiate.
 */
void IonEquationTerm::SetJacobianBlock(
    const len_t derivId, const len_t uqtyId, FVM::Matrix *jac
) {
    const len_t nr = grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 < Zion; Z0++, idx++)
        this->SetCSJacobianBlock(derivId, uqtyId, jac, iIon, Z0, idx*nr);
}

/**
 * Set the elements of the linear operator matrix
 * corresponding to this term.
 * 
 * mat: Linear operator matrix.
 * rhs: Right-hand-side vector.
 */
void IonEquationTerm::SetMatrixElements(
    FVM::Matrix *mat, real_t *rhs
) {
    const len_t nr = grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 < Zion; Z0++, idx++)
        this->SetCSMatrixElements(mat, rhs, iIon, Z0, idx*nr);
}

/**
 * Set the function vector corresponding to this
 * equation term.
 *
 * vec: Function vector.
 * x:   Unknown quantity vector.
 */
void IonEquationTerm::SetVectorElements(
    real_t *vec, const real_t *x
) {
    const len_t nr = grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 < Zion; Z0++, idx++)
        this->SetCSVectorElements(vec, x, iIon, Z0, idx*nr);
}

