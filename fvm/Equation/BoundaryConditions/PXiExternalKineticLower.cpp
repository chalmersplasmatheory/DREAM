/**
 * Reference implementation of the kinetic-kinetic boundary condition on the lower
 * momentum grid.
 */

#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticLower.hpp"


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
PXiExternalKineticLower::PXiExternalKineticLower(
    DREAM::FVM::Grid *lowerGrid, DREAM::FVM::Grid *upperGrid, const DREAM::FVM::Operator *op,
    const len_t id_f_low, const len_t id_f_upp
) : BoundaryCondition(lowerGrid), lowerGrid(lowerGrid), upperGrid(upperGrid), oprtr(op), id_f_low(id_f_low), id_f_upp(id_f_upp) { }

/**
 * Destructor.
 */
PXiExternalKineticLower::~PXiExternalKineticLower() {
}

/**
 * Rebuild the coefficients of this boundary condition.
 */
bool PXiExternalKineticLower::Rebuild(const real_t, DREAM::FVM::UnknownQuantityHandler *unknowns) {
    this->fLow = unknowns->GetUnknownData(id_f_low);
    this->fUpp = unknowns->GetUnknownData(id_f_upp);

    return true;
}

/**
 * Add elements to the Jacobian.
 */
void PXiExternalKineticLower::AddToJacobianBlock(
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
void PXiExternalKineticLower::AddToMatrixElements(
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
void PXiExternalKineticLower::AddToVectorElements(
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
void PXiExternalKineticLower::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> fLow,
    std::function<void(const len_t, const len_t, const real_t)> fUpp
) {
    const len_t NR = this->lowerGrid->GetNr();

    len_t loffset = 0, uoffset = 0;
    for (len_t ir = 0; ir < NR; ir++) {
        const DREAM::FVM::MomentumGrid *mgLower = this->lowerGrid->GetMomentumGrid(ir);
        const DREAM::FVM::MomentumGrid *mgUpper = this->upperGrid->GetMomentumGrid(ir);

        const len_t
            np      = mgLower->GetNp1(),
            nxi     = mgLower->GetNp2(),
            unp     = mgUpper->GetNp1(),
            unxi    = mgUpper->GetNp2();

        const real_t
            *p      = mgLower->GetP1(),
            *up     = mgUpper->GetP1(),
            //*p_f    = mgLower->GetP1_f(),
            //*xi     = mgLower->GetP2(),
            *xi_f   = mgLower->GetP2_f(),
            *Vp_f1  = this->lowerGrid->GetVp_f1(ir),
            *Vp     = this->lowerGrid->GetVp(ir),
            *dp     = mgLower->GetDp1(),
            *dxi    = mgLower->GetDp2(),
            *udxi   = mgUpper->GetDp2(),
            *uxi_f  = mgUpper->GetP2_f(),
            *udp    = mgUpper->GetDp1(),
            *uVp_f1 = this->upperGrid->GetVp_f1(ir);

        // Iterate over xi
        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = loffset + np*j + np-1;

            const real_t D11    = oprtr->GetDiffusionCoeff11(ir)[(np+1)*j + np];
            const real_t F1     = oprtr->GetAdvectionCoeff1(ir)[(np+1)*j + np];
            //const real_t *delta = oprtr->GetInterpolationCoeff1(ir, np, j);
            real_t delta = (F1 > 0 ? 1 : 0);

            real_t aS1 = F1*Vp_f1[j*(np+1) + np] / (Vp[j*np+np-1]*dp[np-1]);
            real_t dS1 = D11*Vp_f1[j*(np+1) + np] / (Vp[j*np+np-1]*dp[np-1]*(up[0]-p[np-1]));
            
            // Advection
            fLow(idx, idx, delta * aS1);

            // Diffusion
            fLow(idx, idx, +dS1);

            // Sum over xi cells
            for (len_t J = 0; J  < unxi; J++) {
                if (!(xi_f[j] <= uxi_f[J+1] && uxi_f[J] <= xi_f[j+1]))
                    continue;

                real_t dxiBar = std::min(xi_f[j+1], uxi_f[J+1]) - std::max(xi_f[j], uxi_f[J]);
                if (dxiBar == 0)
                    continue;

                real_t geom = (dxiBar*uVp_f1[J*(unp+1)]) / (dxi[j]*Vp_f1[j*(np+1)+np]);

                // Advection
                fUpp(idx, uoffset+unp*J, -(1-delta)*aS1*geom);

                // Diffusion
                fUpp(idx, uoffset+unp*J, -dS1*geom);
            }
        }

        loffset += np*nxi;
        uoffset += unp*unxi;
    }
}

