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
 */
void PXiAdvectionDiffusionBoundaryCondition::AddPartialJacobianContributions(
    const len_t, const len_t derivId, Matrix *jac, const real_t *x
) {
    // Iterate over all advection operators...
    const len_t nr = this->grid->GetNr();
    const AdvectionDiffusionTerm *adt = oprtr->GetAdvectionDiffusion();
    for (AdvectionTerm *at : adt->GetAdvectionTerms()) {
        len_t nMultiples;
        if (!at->HasJacobianContribution(derivId, &nMultiples))
            continue;

        at->SetPartialAdvectionTerm(derivId, nMultiples);

        for (len_t n = 0; n < nMultiples; n++)
            SetPartialJacobianContribution(
                n, jac, x,
                at->GetAdvectionDiffCoeff1()+n*nr,
                at->GetAdvectionDiffCoeff2()+n*nr
            );
    }

    // Iterate over all diffusion operators...
    for (DiffusionTerm *dt : adt->GetDiffusionTerms()) {
        len_t nMultiples;
        if (!dt->HasJacobianContribution(derivId, &nMultiples))
            continue;

        dt->SetPartialDiffusionTerm(derivId, nMultiples);

        for (len_t n = 0; n < nMultiples; n++)
            SetPartialJacobianContribution(
                n, jac, x, nullptr, nullptr,
                dt->GetDiffusionDiffCoeff11()+n*nr,
                dt->GetDiffusionDiffCoeff12()+n*nr,
                dt->GetDiffusionDiffCoeff21()+n*nr,
                dt->GetDiffusionDiffCoeff22()+n*nr
            );
    }
}

/**
 * Sets the elements of a block in the jacobian matrix corresponding
 * to a derivative with respect to a *fluid* quantity.
 *
 * n:    Index of multiple to set.
 * jac:  Jacobian matrix to set elements in.
 * x:    Unknown quantity vector (i.e. f).
 * df1:  Derivative of p1-advection coefficient w.r.t. the kinetic quantity
 * dd11: Derivative of p1p1-diffusion coefficient w.r.t. the kinetic quantity
 * dd11: Derivative of p1p2-diffusion coefficient w.r.t. the kinetic quantity
 */
void PXiAdvectionDiffusionBoundaryCondition::SetPartialJacobianContribution(
    const len_t n, Matrix *jac, const real_t *x,
    const real_t *const* df1, const real_t *const*,
    const real_t *const* dd11, const real_t *const* dd12,
    const real_t *const*, const real_t *const*
) {
    ResetJacobianColumn();
    AddToVectorElements_c(jacobianColumn, x, df1, dd11, dd12);

    const len_t nr = this->grid->GetNr();
    for (len_t ir = 0, offset = 0; ir < nr; ir++) {
        // For a fluid grid, these are n1=n2=1
        const len_t n1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t n2 = this->grid->GetMomentumGrid(ir)->GetNp2();

        for (len_t j = 0; j < n2; j++)
            for (len_t i = 0; i < n1; i++)
                jac->SetElement(offset + n1*j + i, n*nr+ir, jacobianColumn[offset + n1*j + i]);
        
        offset += n1*n2;
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
    const real_t *const* Ap  = oprtr->GetAdvectionCoeff1();
    const real_t *const* Ax  = oprtr->GetAdvectionCoeff2();
    const real_t *const* Dpp = oprtr->GetDiffusionCoeff11();
    const real_t *const* Dpx = oprtr->GetDiffusionCoeff12();
    const real_t *const* Dxp = oprtr->GetDiffusionCoeff21();
    const real_t *const* Dxx = oprtr->GetDiffusionCoeff22();

    this->AddToVectorElements_c(vec, f, Ap, Ax, Dpp, Dpx, Dxp, Dxx);
}

