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
    
    SetName("TransportBC");

    if (type == TRANSPORT_BC_DF_CONST && grid->GetNr() == 1)
        throw DREAMException(
            "Transport boundary condition: The 'DF_CONST' boundary condition "
            "can only be applied to radial grids with nr > 1."
        );
    
    this->jacobianColumn = new real_t[grid->GetNCells()];
}

/**
 * Destructor.
 */
template<typename T>
DREAM::TransportBC<T>::~TransportBC() {
    delete [] this->jacobianColumn;
}

/**
 * Called when the grid object is rebuilt.
 */
template<typename T>
bool DREAM::TransportBC<T>::GridRebuilt() {
    delete [] this->jacobianColumn;

    this->jacobianColumn = new real_t[this->grid->GetNCells()];

    return true;
}

/**
 * Add elements to the given Jacobian.
 */
template<typename T>
bool DREAM::TransportBC<T>::AddToJacobianBlock(
    const len_t uqtyId, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t *x
) {
    bool contributes = (derivId == uqtyId);
    // Derivative w.r.t. the quantity which the operator is applied to?
    if (derivId == uqtyId)
        this->AddToMatrixElements(jac, nullptr);

    const len_t nr = this->grid->GetNr();

    // Add off-diagonal derivatives
    len_t nMultiples;
    if (!transportOperator->HasJacobianContribution(derivId, &nMultiples))
        return contributes;

    SetPartialTerm(derivId, nMultiples);

    for (len_t n = 0; n < nMultiples; n++) {
        SetPartialJacobianContribution(
            0, n, jac, x, JACOBIAN_SET_CENTER,
            GetDiffCoefficient()+n*(nr+1)
        );

        // Add terms corresponding to radial interpolation
        // (for coefficients evaluated on the flux grid, which
        // depend on quantities only known on the radial
        // distribution grid)
        SetPartialJacobianContribution(
            -1, n, jac, x, JACOBIAN_SET_LOWER,
            GetDiffCoefficient()+n*(nr+1)
        );

        SetPartialJacobianContribution(
            +1, n, jac, x, JACOBIAN_SET_UPPER,
            GetDiffCoefficient()+n*(nr+1)
        );
    }

    return contributes;
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
    const len_t 
        nr = this->grid->GetNr(),
        ir = nr-1;

    this->__SetElements([&mat](const len_t I, const len_t J, const real_t V) {
        mat->SetElement(I, J, V);
    }, GetCoefficient(ir+1), NO_JACOBIAN);
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
    AddToVectorElements_c(vec, f, GetCoefficient(), NO_JACOBIAN);
}

template<typename T>
void DREAM::TransportBC<T>::AddToVectorElements_c(
    real_t *vec, const real_t *f, const real_t *const* coeff, jacobian_interp_mode set_mode
) {
    const len_t 
        nr = this->grid->GetNr(),
        ir = nr-1;

    this->__SetElements([&vec,&f](const len_t I, const len_t J, const real_t V) {
        vec[I] += V*f[J];
    }, coeff[ir+1], set_mode);
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
    std::function<void(const len_t, const len_t, const real_t)> f,
    jacobian_interp_mode set_mode
) {
    const len_t 
        nr = this->grid->GetNr(),
        ir = nr-1;

    __SetElements(f, GetCoefficient(ir+1));
}

template<typename T>
void DREAM::TransportBC<T>::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> f,
    const real_t *coeff, jacobian_interp_mode set_mode
) {
    const len_t 
        nr = this->grid->GetNr(),
        ir = nr-1, 
        np1 = grid->GetNp1(ir),
        np2 = grid->GetNp2(ir),
        offset = grid->GetNCells() - np1*np2;
        
    const real_t
        *Vp_fr = this->grid->GetVp_fr(ir+1),
        *Vp    = this->grid->GetVp(ir),
        dr     = this->grid->GetRadialGrid()->GetDr(ir),
        dr1    = this->grid->GetRadialGrid()->GetDr(ir-1),  // out-of-bounds checked for in constructor
        *deltaRadialFlux = this->transportOperator->GetRadialJacobianInterpolationCoeffs();

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

        // Interpolation in radius
        if (set_mode == JACOBIAN_SET_LOWER)
            v = 0;
        else if (set_mode == JACOBIAN_SET_CENTER)
            v *= 1 - deltaRadialFlux[ir];
        else if (set_mode == JACOBIAN_SET_UPPER)
            v *= deltaRadialFlux[ir];

        // Select appropriate form for the boundary condition
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

/**
 * Sets the elements of a block in the jacobian matrix corresponding
 * to a derivative with respect to a *fluid* quantity.
 *
 * diagOffs: A non-zero value indicates that a sub- or super-diagonal should
 *           be set, rather than the main diagonal.
 * n:        Index of multiple to set.
 * jac:      Jacobian matrix to set elements in.
 * x:        Unknown quantity vector.
 * set_mode: Indicates how to set the elements (for when radial interpolation is
 *           necessary).
 * dfr:      Derivative of r-advection coefficient w.r.t. the unknown quantity.
 * ddrr:     Derivative of rr-diffusion coefficient w.r.t. the unknown quantity.
 */
template<typename T>
void DREAM::TransportBC<T>::SetPartialJacobianContribution(
    const int_t diagOffs, const len_t n, DREAM::FVM::Matrix *jac, const real_t *x,
    jacobian_interp_mode set_mode, const real_t *const* coeff
) {
    ResetJacobianColumn();
    AddToVectorElements_c(jacobianColumn, x, coeff, set_mode);

    const len_t nr = this->grid->GetNr();
    for (len_t ir = 0, offset = 0; ir < nr; ir++) {
        if ((ir==0 && diagOffs==-1) || (ir+diagOffs >= nr))
            continue;

        const len_t N = this->grid->GetMomentumGrid(ir)->GetNCells();
        len_t col = n*nr + ir + diagOffs;
        for (len_t i = offset; i < offset+N; i++)
            jac->SetElement(i, col, jacobianColumn[i]);

        offset += N;
    }
}

/**
 * Reset the vector that is used to temporarily store a single
 * column of the jacobian when building partial contributions
 * to the jacobian matrix.
 */
template<typename T>
void DREAM::TransportBC<T>::ResetJacobianColumn() {
    const len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        jacobianColumn[i] = 0;
}

