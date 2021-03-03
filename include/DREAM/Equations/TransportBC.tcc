/**
 * Implementation of a boundary condition for an advection/diffusion term
 * in the prescribed transport model of DREAM.
 */

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/TransportBC.hpp"


/**
 * Constructor.
 */
template<typename T>
DREAM::TransportBC<T>::TransportBC(
    DREAM::FVM::Grid *grid, T *tt,
    enum bctype type
) : FVM::BC::BoundaryCondition(grid), transportOperator(tt), type(type) {

    if (type == TRANSPORT_BC_DF_CONST && grid->GetNr() == 1)
        throw DREAMException(
            "Transport boundary condition: The 'DF_CONST' boundary condition "
            "can only be applied to radial grids with nr > 1."
        );
}


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
        dr     = this->grid->GetRadialGrid()->GetDr(ir),
        dr1    = this->grid->GetRadialGrid()->GetDr(ir-1);  // out-of-bounds checked for in constructor

    real_t dr_f;
    if (ir == 0)
        dr_f = 1;
    else
        dr_f = this->grid->GetRadialGrid()->GetDr_f(ir-1);

    // Iterate over every momentum cell...
    const real_t Nm = np1*np2;
    for (len_t idx = 0; idx < Nm; idx++) {

        // Flux (without advection/diffusion coefficient)
        real_t S_wo_coeff =
            Vp_fr[idx] / (Vp[idx] * dr);

        real_t v = __GetSingleElement(coeff[idx], S_wo_coeff, dr_f);

        switch (this->type) {
            case TRANSPORT_BC_F0:
                f(offset+idx, offset+idx, v);
                break;

            case TRANSPORT_BC_DF_CONST: {
                f(offset+idx, offset+idx, v);

                // Set T_{N+1} = T_N + dr_N * (T_N - T_{N-1}) / dr_{N-1}
                //               = (1+delta)*T_N - delta*T_{N-1}
                real_t delta = dr / dr1;
                f(offset+idx, offset+idx, -v*(1+delta));
                f(offset+idx, offset-Nm+idx, v*delta);
            } break;

            default:
                throw DREAMException("Unrecognized transport boundary condition specified: %d.", this->type);
        }
    }
}
