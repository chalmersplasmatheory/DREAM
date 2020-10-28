/**
 * Implementation of a boundary condition for an advection/diffusion term
 * in the prescribed transport model of DREAM.
 */

#include "DREAM/Equations/TransportBC.hpp"


/**
 * Constructor.
 */
template<typename T>
DREAM::TransportBC<T>::TransportBC(
    DREAM::FVM::Grid *grid, T *tt
) : FVM::BC::BoundaryCondition(grid), transportOperator(tt) { }


/**
 * Add elements to the given Jacobian.
 */
template<typename T>
void DREAM::TransportBC<T>::AddToJacobianBlock(
    const len_t uqtyId, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t*
) {
    if (derivId == uqtyId)
        this->AddToMatrixElements(jac, nullptr);

    // TODO handle derivatives of coefficients
}

/**
 * Add flux to linearized operator matrix.
 *
 * mat: Matrix to add boundary condition to.
 * rhs: Right-hand-side vector (not used).
 */
template<typename T>
void DREAM::TransportBC<T>::AddToMatrixElements(
    DREAM::FVM::Matrix *mat, real_t*
) {
    this->__SetElements([&mat](const len_t I, const len_t J, const real_t V) {
        mat->SetElement(I, J, V);
    });
}

/**
 * Add flux to function vector.
 *
 * vec: Function vector to add flux to.
 * x:   Vector of values for unknown quantity to which this operator is applied.
 */
template<typename T>
void DREAM::TransportBC<T>::AddToVectorElements(
    real_t *vec, const real_t *f
) {
    this->__SetElements([&vec,&f](const len_t I, const len_t J, const real_t V) {
        vec[I] += V*f[J];
    });
}

/**
 * PRIVATE
 * Internal routine used for setting matrix/vector elements.
 *
 * f(I,J,V): Function for setting matrix/vector elements. I denotes the
 *           index of the unknown quantity to set (matrix row), J denotes
 *           the inedx of the distribution function to evaluate, and V is
 *           a scalar value to weight the distribution function with.
 */
template<typename T>
void DREAM::TransportBC<T>::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> f
) {
    const len_t 
        nr = this->grid->GetNr(),
        ir = nr-1, 
        np1 = grid->GetNp1(ir),
        np2 = grid->GetNp2(ir),
        offset = grid->GetNCells() - np1*np2;
        
    const real_t *coeff = this->GetCoefficient(ir+1);

    const real_t
        *Vp_fr = this->grid->GetVp_fr(ir+1),
        *Vp    = this->grid->GetVp(ir),
        dr     = this->grid->GetRadialGrid()->GetDr(ir);

    real_t dr_f;
    if (ir == 0)
        dr_f = 1;
    else
        dr_f = this->grid->GetRadialGrid()->GetDr_f(ir-1);

    // Iterate over every momentum cell...
    for (len_t idx = 0; idx < np1*np2; idx++) {

        // Flux (without advection/diffusion coefficient)
        real_t S_wo_coeff =
            Vp_fr[idx] / (Vp[idx] * dr);

        real_t v = __GetSingleElement(coeff[idx], S_wo_coeff, dr_f);
        f(offset+idx, offset+idx, v);
    }
}
