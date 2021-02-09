/**
 * This class implements some general routines for building jacobian
 * matrices for advection-diffusion boundary conditions on p/xi grids.
 * The purpose of this class is for the more specialized boundary conditions
 * to be derived from it.
 */

#include "FVM/Equation/BoundaryConditions/PXiAdvectionDiffusionBoundaryCondition.hpp"


using namespace DREAM::FVM::BC;


/**
 * Constructor.
 */
PXiAdvectionDiffusionBoundaryCondition::
    PXiAdvectionDiffusionBoundaryCondition(Grid *g, const Operator *o)
        : BoundaryCondition(g), oprtr(o) {

    this->jacobianColumn = new real_t[g->GetNCells()];
}

/**
 * Destructor.
 */
PXiAdvectionDiffusionBoundaryCondition::
    ~PXiAdvectionDiffusionBoundaryCondition() {

    delete [] this->jacobianColumn;
}

/**
 * Method which is called whenever the parent grid has
 * been rebuilt.
 */
bool PXiAdvectionDiffusionBoundaryCondition::GridRebuilt() {
    delete [] this->jacobianColumn;

    this->jacobianColumn = new real_t[this->grid->GetNCells()];

    return true;
}

/**
 * Add elements to the jacobian matrix for this boundary condition.
 * NOTE: This method DOES NOT add the usual diagonal contribution resulting
 * differentiating the quantity 'qtyId' with respect to itself. It only
 * handles the terms which contain derivatives of the advection and
 * diffusion coefficients!
 *
 * qtyId:   ID of the unknown quantity which should be differentiated.
 * derivId: ID of the unknown quantity w.r.t. which the operator should be differentiated.
 * jac:     Jacobian matrix to add elements into.
 * x:       Current value of the unknown quantity 'qtyId'.
 *
 * hasRadialContributions: When 'true', this flag indicates that this boundary
 *                         conditions has coefficients which are to be evaluated
 *                         on the radial _flux_ grid, and should thus be
 *                         interpolated across several points on the distribution
 *                         grid.
 */
void PXiAdvectionDiffusionBoundaryCondition::AddPartialJacobianContributions(
    const len_t, const len_t derivId, Matrix *jac, const real_t *x,
    bool hasRadialContributions
) {
    // Iterate over all advection operators...
    const len_t nr = this->grid->GetNr();
    const AdvectionDiffusionTerm *adt = oprtr->GetAdvectionDiffusion();
    for (AdvectionTerm *at : adt->GetAdvectionTerms()) {
        len_t nMultiples;
        if (!at->HasJacobianContribution(derivId, &nMultiples))
            continue;

        at->SetPartialAdvectionTerm(derivId, nMultiples);

        for (len_t n = 0; n < nMultiples; n++) {
            SetPartialJacobianContribution(
                0, n, jac, x, JACOBIAN_SET_CENTER,
                at->GetAdvectionDiffCoeffR()+n*(nr+1),
                at->GetAdvectionDiffCoeff1()+n*nr,
                at->GetAdvectionDiffCoeff2()+n*nr,
                nullptr, nullptr, nullptr,
                nullptr, nullptr
            );

            // Add terms corresponding to radial interpolation
            // (for coefficients evaluated on the flux grid, which
            // depend on quantities only known on the radial
            // distribution grid)
            if (hasRadialContributions) {
                SetPartialJacobianContribution(
                    -1, n, jac, x, JACOBIAN_SET_LOWER,
                    at->GetAdvectionDiffCoeffR()+n*(nr+1),
                    at->GetAdvectionDiffCoeff1()+n*nr,
                    at->GetAdvectionDiffCoeff2()+n*nr,
                    nullptr, nullptr, nullptr,
                    nullptr, nullptr
                );

                SetPartialJacobianContribution(
                    +1, n, jac, x, JACOBIAN_SET_UPPER,
                    at->GetAdvectionDiffCoeffR()+n*(nr+1),
                    at->GetAdvectionDiffCoeff1()+n*nr,
                    at->GetAdvectionDiffCoeff2()+n*nr,
                    nullptr, nullptr, nullptr,
                    nullptr, nullptr
                );
            }
        }
    }

    // Iterate over all diffusion operators...
    for (DiffusionTerm *dt : adt->GetDiffusionTerms()) {
        len_t nMultiples;
        if (!dt->HasJacobianContribution(derivId, &nMultiples))
            continue;

        dt->SetPartialDiffusionTerm(derivId, nMultiples);

        for (len_t n = 0; n < nMultiples; n++) {
            SetPartialJacobianContribution(
                0, n, jac, x, JACOBIAN_SET_CENTER,
                nullptr, nullptr, nullptr,
                dt->GetDiffusionDiffCoeffRR()+n*(nr+1),
                dt->GetDiffusionDiffCoeff11()+n*nr,
                dt->GetDiffusionDiffCoeff12()+n*nr,
                dt->GetDiffusionDiffCoeff21()+n*nr,
                dt->GetDiffusionDiffCoeff22()+n*nr
            );

            if (hasRadialContributions) {
                SetPartialJacobianContribution(
                    -1, n, jac, x, JACOBIAN_SET_LOWER,
                    nullptr, nullptr, nullptr,
                    dt->GetDiffusionDiffCoeffRR()+n*(nr+1),
                    dt->GetDiffusionDiffCoeff11()+n*nr,
                    dt->GetDiffusionDiffCoeff12()+n*nr,
                    dt->GetDiffusionDiffCoeff21()+n*nr,
                    dt->GetDiffusionDiffCoeff22()+n*nr
                );

                SetPartialJacobianContribution(
                    +1, n, jac, x, JACOBIAN_SET_UPPER,
                    nullptr, nullptr, nullptr,
                    dt->GetDiffusionDiffCoeffRR()+n*(nr+1),
                    dt->GetDiffusionDiffCoeff11()+n*nr,
                    dt->GetDiffusionDiffCoeff12()+n*nr,
                    dt->GetDiffusionDiffCoeff21()+n*nr,
                    dt->GetDiffusionDiffCoeff22()+n*nr
                );
            }
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
 * x:        Unknown quantity vector (i.e. f).
 * set_mode: Indicates how to set the elements (for when radial interpolation is
 *           necessary).
 * df1:      Derivative of p1-advection coefficient w.r.t. the kinetic quantity
 * dd11:     Derivative of p1p1-diffusion coefficient w.r.t. the kinetic quantity
 * dd11:     Derivative of p1p2-diffusion coefficient w.r.t. the kinetic quantity
 */
void PXiAdvectionDiffusionBoundaryCondition::SetPartialJacobianContribution(
    const int_t diagOffs, const len_t n, Matrix *jac, const real_t *x,
    jacobian_interp_mode set_mode,
    const real_t *const* dfr, const real_t *const* df1, const real_t *const* df2,
    const real_t *const* ddrr, const real_t *const* dd11, const real_t *const* dd12,
    const real_t *const* dd21, const real_t *const* dd22
) {
    ResetJacobianColumn();
    AddToVectorElements_c(jacobianColumn, x, dfr, df1, df2, ddrr, dd11, dd12, dd21, dd22, set_mode);

    const len_t nr = this->grid->GetNr();
    for (len_t ir = 0, offset = 0; ir < nr; ir++) {
        if ((ir==0 && diagOffs==-1) || (ir+diagOffs >= nr))
            continue;

        // For a fluid grid, N=1
        const len_t N = this->grid->GetMomentumGrid(ir)->GetNCells();
        len_t col = n*nr+ir+diagOffs;
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
void PXiAdvectionDiffusionBoundaryCondition::ResetJacobianColumn() {
    const len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        jacobianColumn[i] = 0;
}

/**
 * Add flux to function vector.
 */
void PXiAdvectionDiffusionBoundaryCondition::AddToVectorElements(
    real_t *vec, const real_t *f
) {
    const real_t *const* Ar  = oprtr->GetAdvectionCoeffR();
    const real_t *const* Ap  = oprtr->GetAdvectionCoeff1();
    const real_t *const* Ax  = oprtr->GetAdvectionCoeff2();
    const real_t *const* Drr = oprtr->GetDiffusionCoeffRR();
    const real_t *const* Dpp = oprtr->GetDiffusionCoeff11();
    const real_t *const* Dpx = oprtr->GetDiffusionCoeff12();
    const real_t *const* Dxp = oprtr->GetDiffusionCoeff21();
    const real_t *const* Dxx = oprtr->GetDiffusionCoeff22();

    this->AddToVectorElements_c(vec, f, Ar, Ap, Ax, Drr, Dpp, Dpx, Dxp, Dxx);
}
