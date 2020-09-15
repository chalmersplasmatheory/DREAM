/**
 * Reference implementation of the kinetic-kinetic boundary condition on the upper
 * momentum grid.
 */

#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticUpper.hpp"


using namespace DREAM::FVM::BC;


/**
 * Constructor.
 *
 * lowerGrid: Grid on which the lower distribution function lives.
 * upperGrid: Grid on which the upper distribution function lives.
 * op:        Advection-diffusion operator applied to the lower distribution function.
 * id_f_low:  ID of the lower distribution function.
 * id_f_upp:  ID of the upper distribution function.
 */
PXiExternalKineticUpper::PXiExternalKineticUpper(
    DREAM::FVM::Grid *lowerGrid, DREAM::FVM::Grid *upperGrid, const DREAM::FVM::Operator *op,
    const len_t id_f_low, const len_t id_f_upp
) : BoundaryCondition(upperGrid), lowerGrid(lowerGrid), upperGrid(upperGrid), oprtr(op), id_f_low(id_f_low), id_f_upp(id_f_upp) { }

/**
 * Destructor.
 */
PXiExternalKineticUpper::~PXiExternalKineticUpper() {
}

/**
 * Rebuild the coefficients of this boundary condition.
 */
bool PXiExternalKineticUpper::Rebuild(const real_t, DREAM::FVM::UnknownQuantityHandler *unknowns) {
    this->fLow = unknowns->GetUnknownData(id_f_low);
    this->fUpp = unknowns->GetUnknownData(id_f_upp);

    return true;
}

/**
 * Add elements to the Jacobian.
 */
void PXiExternalKineticUpper::AddToJacobianBlock(
    const len_t uqtyId, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t*
) {
    if (derivId == uqtyId)
        this->AddToMatrixElements(jac, nullptr);

    // TODO handle derivatives of coefficients...
}

/**
 * Add flux to linearized operator matrix.
 * 
 * mat: Matrix to add boundary conditions to.
 * rhs: Right-hand-side vector (not used).
 */
void PXiExternalKineticUpper::AddToMatrixElements(
    Matrix *mat, real_t*
) {
    // Reset matrix offsets...
    const len_t rowOffset = mat->GetRowOffset();
    const len_t colOffset = mat->GetColOffset();
    mat->SetOffset(rowOffset, 0);

    // Get matrix offsets of f_low and f_upp
    BlockMatrix *bm = static_cast<BlockMatrix*>(mat);
    len_t lowOffset = bm->GetOffsetById(id_f_low);
    len_t uppOffset = bm->GetOffsetById(id_f_upp);

    __SetElements(
        [mat,lowOffset](const len_t I, const len_t J, const real_t V) { mat->SetElement(I, lowOffset+J, V); },
        [mat,uppOffset](const len_t I, const len_t J, const real_t V) { mat->SetElement(I, uppOffset+J, V); }
    );

    // Restore matrix offset
    mat->SetOffset(rowOffset, colOffset);
}

/**
 * Add flux to function vector.
 */
void PXiExternalKineticUpper::AddToVectorElements(
    real_t *vec, const real_t*
) {
    const real_t *fLow = this->fLow;
    const real_t *fUpp = this->fUpp;

    __SetElements(
        [vec,fLow](const len_t I, const len_t J, const real_t V) { vec[I] += V*fLow[J]; },
        [vec,fUpp](const len_t I, const len_t J, const real_t V) { vec[I] += V*fUpp[J]; }
    );
}

/**
 * Set matrix or vector elements.
 */
void PXiExternalKineticUpper::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> fLow,
    std::function<void(const len_t, const len_t, const real_t)> fUpp
) {
    const len_t NR = this->upperGrid->GetNr();

    len_t loffset = 0, uoffset = 0;
    for (len_t ir = 0; ir < NR; ir++) {
        const DREAM::FVM::MomentumGrid *mgLower = this->lowerGrid->GetMomentumGrid(ir);
        const DREAM::FVM::MomentumGrid *mgUpper = this->upperGrid->GetMomentumGrid(ir);

        const len_t
            np   = mgUpper->GetNp1(),
            nxi  = mgUpper->GetNp2(),
            lnp  = mgLower->GetNp1(),
            lnxi = mgLower->GetNp2();

        const real_t
            *p      = mgUpper->GetP1(),
            *lp     = mgLower->GetP1(),
            *xi_f   = mgUpper->GetP2_f(),
            *lxi_f  = mgLower->GetP2_f(),
            *Vp_f1  = this->upperGrid->GetVp_f1(ir),
            *Vp     = this->upperGrid->GetVp(ir),
            *lVp_f  = this->lowerGrid->GetVp_f1(ir),
            *dp     = mgUpper->GetDp1(),
            *dxi    = mgUpper->GetDp2(),
            *ldp    = mgLower->GetDp1(),
            *ldxi   = mgLower->GetDp2();

        // Iterate over xi
        for (len_t J = 0; J < nxi; J++) {
            const len_t idx = uoffset + np*J + 0;

            for (len_t j = 0; j < lnxi; j++) {
                const len_t lidx = loffset + lnp*j + np-1;

                if (!(xi_f[J] <= lxi_f[j+1] && lxi_f[j] <= xi_f[J+1]))
                    continue;

                real_t dxiBar = std::min(lxi_f[j+1], xi_f[J+1]) - std::max(lxi_f[j], xi_f[J]);
                if (dxiBar == 0)
                    continue;

                //const real_t *delta = oprtr->GetInterpolationCoeff1(ir, np, j);
                const real_t D11 = oprtr->GetDiffusionCoeff11(ir)[(lnp+1)*j + lnp];
                const real_t F1  = oprtr->GetAdvectionCoeff1(ir)[(lnp+1)*j + lnp];
                const real_t delta = (F1 >= 0 ? 1 : 0);

                real_t aS1 = F1*Vp_f1[J*(np+1)] / (Vp[J*np]*dp[0]);
                real_t dS1 = D11*Vp_f1[J*(np+1)] / (Vp[J*np]*dp[0]*(p[0]-lp[lnp-1]));

                real_t geom = (dxiBar*lVp_f[j*(lnp+1)+lnp]) / (dxi[J]*Vp_f1[J*(np+1)]);

                // Advection
                fLow(idx, lidx, -delta*aS1*geom);
                fUpp(idx, idx,   (1-delta)*aS1*dxiBar/ldxi[j]);

                // Diffusion
                fLow(idx, lidx, -dS1*geom);
                fUpp(idx, idx,  +dS1*dxiBar/ldxi[j]);
            }
        }

        loffset += lnp*lnxi;
        uoffset += np*nxi;
    }
}

