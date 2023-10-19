/**
 * Implementation of a general ion boundary condition, which applies to
 * all ions and ion charge states in the equation system.
 */

#include "DREAM/Equations/Fluid/IonBoundaryCondition.hpp"


/**
 * Constructor.
 */
template<class T>
template<typename ... Args>
IonBoundaryCondition<T>::IonBoundaryCondition(
    FVM::Grid *g, IonHandler *ihdl, const len_t iIon,
    Args&& ... args
) : T(g, std::forward<Args>(args) ...), ions(ihdl), iIon(iIon) {

    this->Zion = this->ions->GetZ(iIon);
}

/**
 * Destructor.
 */
template<class T>
IonBoundaryCondition<T>::~IonBoundaryCondition() {}


/**
 * Add elements to the specified block in the given Jacobian matrix.
 *
 * derivId: ID of the unknown quantity with respect to
 *          which differentiation should be done.
 * uqtyId:  ID of the unknown quantity to differentiate.
 */
template<class T>
bool IonBoundaryCondition<T>::AddToJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    FVM::Matrix *jac, const real_t *x
) {
    bool contributes = false;
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        contributes |=
            this->AddCSToJacobianBlock(uqtyId, derivId, jac, x, iIon, Z0, idx*nr);

    return contributes;
}

/**
 * Set elements of the specified block in the given Jacobian matrix.
 *
 * derivId: ID of the unknown quantity with respect to
 *          which differentiation should be done.
 * uqtyId:  ID of the unknown quantity to differentiate.
 */
template<class T>
bool IonBoundaryCondition<T>::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    FVM::Matrix *jac, const real_t *x
) {
    bool contributes = false;
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        contributes |=
            this->SetCSJacobianBlock(uqtyId, derivId, jac, x, iIon, Z0, idx*nr);

    return contributes;
}

/**
 * Add to the elements of the linear operator matrix corresponding
 * to this term.
 *
 * mat: Linear operator matrix.
 * rhs: Right-hand side vector.
 */
template<class T>
void IonBoundaryCondition<T>::AddToMatrixElements(
    FVM::Matrix *mat, real_t *rhs
) {
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        this->AddCSToMatrixElements(mat, rhs, iIon, Z0, idx*nr);
}

/**
 * Set the elements of the linear operator matrix corresponding
 * to this term.
 *
 * mat: Linear operator matrix.
 * rhs: Right-hand side vector.
 */
template<class T>
void IonBoundaryCondition<T>::SetMatrixElements(
    FVM::Matrix *mat, real_t *rhs
) {
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        this->SetCSMatrixElements(mat, rhs, iIon, Z0, idx*nr);
}

/**
 * Add to the function vector corresponding to this
 * equation term.
 *
 * vec: Function vector.
 * x:   Unknown quantity vector.
 */
template<class T>
void IonBoundaryCondition<T>::AddToVectorElements(
    real_t *vec, const real_t *x
) {
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        this->AddCSToVectorElements(vec, x, iIon, Z0, idx*nr);
}

/**
 * Set the function vector corresponding to this
 * equation term.
 *
 * vec: Function vector.
 * x:   Unknown quantity vector.
 */
template<class T>
void IonBoundaryCondition<T>::SetVectorElements(
    real_t *vec, const real_t *x
) {
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for (len_t Z0 = 0; Z0 <= Zion; Z0++, idx++)
        this->SetCSVectorElements(vec, x, iIon, Z0, idx*nr);
}

