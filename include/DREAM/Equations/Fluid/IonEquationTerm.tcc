/**
 * Implementation of a general ion equation term, which applies to
 * all ions and ion charge states in the equation system.
 */

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"



/**
 * Constructor.
 */
template<class T>
IonEquationTerm<T>::IonEquationTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon)
    : T(g), ions(ihdl), iIon(iIon) {

    this->Zion = this->ions->GetZ(iIon);
}


/**
 * Constructor for MomentQuantity.
 */
template<class T>
IonEquationTerm<T>::IonEquationTerm(
    FVM::Grid *momentGrid, FVM::Grid *fGrid, const len_t momentId, const len_t fId, 
    FVM::UnknownQuantityHandler *u, real_t pThreshold, FVM::MomentQuantity::pThresholdMode pMode,
    IonHandler *ihdl, const len_t iIon
) : T(momentGrid, fGrid, momentId, fId, u, pThreshold, pMode), ions(ihdl), iIon(iIon) {
    this->Zion = this->ions->GetZ(iIon);
}


/**
 * Destructor.
 */
template<class T>
IonEquationTerm<T>::~IonEquationTerm() {}


/**
 * Set the specified block in the given Jacobian matrix.
 *
 * derivId: ID of the unknown quantity with respect to
 *          which differentiation should be done.
 * uqtyId:  ID of the unknown quantity to differentiate.
 */
template<class T>
void IonEquationTerm<T>::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x
) {
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        this->SetCSJacobianBlock(uqtyId, derivId, jac, x, iIon, Z0, idx*nr);
}

/**
 * Set the elements of the linear operator matrix
 * corresponding to this term.
 * 
 * mat: Linear operator matrix.
 * rhs: Right-hand-side vector.
 */
template<class T>
void IonEquationTerm<T>::SetMatrixElements(
    FVM::Matrix *mat, real_t *rhs
) {
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        this->SetCSMatrixElements(mat, rhs, iIon, Z0, idx*nr);
}

/**
 * Set the function vector corresponding to this
 * equation term.
 *
 * vec: Function vector.
 * x:   Unknown quantity vector.
 */
template<class T>
void IonEquationTerm<T>::SetVectorElements(
    real_t *vec, const real_t *x
) {
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        this->SetCSVectorElements(vec, x, iIon, Z0, idx*nr);
}

