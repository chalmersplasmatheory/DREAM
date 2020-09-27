/**
 * Boundary condition representing momentum-space transport across a boundary
 * between two separate grids. The upper flux boundary of the first grid must
 * be equal to the lower flux boundary of the second grid, but the resolution
 * in \xi may be different.
 *
 * This boundary condition can be applied to either of the momentum space grids.
 *      ________________ ________________
 *  ^  |                x                |
 *  |  |     GRID 1     x     GRID 2     |
 * xi  |________________x________________|
 *     p0            p ---->            pmax
 *
 * (the boundary condition is applied on the boundary marked with the crosses)
 */

#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"


using namespace DREAM::FVM::BC;


/**
 * Constructor.
 */
PXiExternalKineticKinetic::PXiExternalKineticKinetic(
    DREAM::FVM::Grid *grid, DREAM::FVM::Grid *lowerGrid, DREAM::FVM::Grid *upperGrid,
    const DREAM::FVM::Operator *eqn, const len_t id_f_low, const len_t id_f_upp,
    enum condition_type ctype
) : BoundaryCondition(grid), lowerGrid(lowerGrid), upperGrid(upperGrid),
    equation(eqn), id_f_low(id_f_low), id_f_upp(id_f_upp), type(ctype) { }

/**
 * Destructor.
 */
PXiExternalKineticKinetic::~PXiExternalKineticKinetic() { }


/**
 * Construct advection/diffusion coefficients to use on the upper grid.
 * These are constructed by averaging the coefficients on the lower grid
 * over the range in \xi covered by each grid cell on the upper grid.
 */
bool PXiExternalKineticKinetic::Rebuild(const real_t, UnknownQuantityHandler *uqh) {
    this->fLow = uqh->GetUnknownData(id_f_low);
    this->fUpp = uqh->GetUnknownData(id_f_upp);
    
    return true;
}

/**
 * Add flux to jacobian block.
 */
void PXiExternalKineticKinetic::AddToJacobianBlock(
    const len_t derivId, const len_t uqtyId, Matrix *jac, const real_t* /*x*/
) {
    if (derivId == uqtyId)
        this->AddToMatrixElements(jac, nullptr);

    // TODO handle derivatives of coefficients
}

/**
 * Add flux to linearized operator matrix.
 * 
 * mat: Matrix to add boundary conditions to.
 * rhs: Right-hand-side vector (not used).
 */
void PXiExternalKineticKinetic::AddToMatrixElements(
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
void PXiExternalKineticKinetic::AddToVectorElements(
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
 * Internal routine for setting matrix/vector elements.
 */
void PXiExternalKineticKinetic::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> fLow,
    std::function<void(const len_t, const len_t, const real_t)> fUpp
) {
    const len_t nr = this->grid->GetNr();
    len_t loffset = 0, uoffset = 0;

    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *lmg = this->lowerGrid->GetMomentumGrid(ir);
        const MomentumGrid *umg = this->upperGrid->GetMomentumGrid(ir);

        const len_t
            lnp = lmg->GetNp1(), lnxi = lmg->GetNp2(),
            unp = umg->GetNp1(), unxi = umg->GetNp2();

        const real_t
            *lp    = lmg->GetP1(),
            *up    = umg->GetP1(),
            *lxi   = lmg->GetP2(),
            *uxi   = umg->GetP2(),
            *lxi_f = lmg->GetP2_f(),
            *uxi_f = umg->GetP2_f(),
            *ldp   = lmg->GetDp1(),
            *udp   = umg->GetDp1(),
            *ldxi  = lmg->GetDp2(),
            *udxi  = umg->GetDp2();

        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        //const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        const real_t
            *lVp   = this->lowerGrid->GetVp(ir),
            *lVp_f = this->lowerGrid->GetVp_f1(ir),
            *uVp   = this->upperGrid->GetVp(ir),
            *uVp_f = this->upperGrid->GetVp_f1(ir),
            *VpVol = this->grid->GetVpVol();

        // j = xi index on lower grid
        // J = xi index on upper grid
        for (len_t j = 0, J = 0; j < lnxi; j++) {

            // Locate correct J (first J such that uxi[J] >= lxi[j])...
            while (uxi[J] < lxi[j] && J < unxi-1)
                J++;

            // Shortcuts to indices...
            len_t
                lidx   = j*lnp+lnp-1,
                lidx_f = j*(lnp+1)+lnp,
                uidx   = J*unp,
                uidx_m = (J-1)*unp,
                uidx_f = J*(unp+1),
                fidx;

            // Set indices for f(r,p,xi) and FVM geometric factor
            // based on which quantity we're building the flux
            // for (for f_hot, f_RE or n_RE)
            real_t Vd;
            if (this->type == TYPE_LOWER) {
                fidx = loffset + lidx;
                Vd   = lVp_f[lidx_f] / (lVp[lidx] * ldp[lnp-1]);
            } else if (this->type == TYPE_UPPER) {
                fidx = uoffset + uidx;
                Vd   = -uVp_f[uidx_f] / (uVp[uidx] * udp[0]);
            } else if (this->type == TYPE_DENSITY) {
                fidx = ir;
                Vd   = -lVp_f[lidx_f] * ldxi[j] / VpVol[ir];
            }

            real_t dxiBar = std::min(lxi_f[j+1], uxi_f[J+1]) - std::max(lxi_f[j], uxi_f[J]);
            real_t fac=1;
            if (this->type == TYPE_LOWER || this->type == TYPE_DENSITY) {
                // XXX This works, but WHY???
                fac = udxi[J] / ldxi[j];
            } else if (this->type == TYPE_UPPER) {
                fac = lVp_f[lidx_f] * dxiBar / (uVp_f[uidx_f] * udxi[J]);
            }

            // Interpolation coefficients...
            real_t delta1 = J==0 ? 1.0 : (lxi[j]-uxi[J-1]) / (uxi[J]-uxi[J-1]);
            real_t delta2 = Ap[lidx_f] > 0 ? 1.0 : 0.0;

            //////////////////////
            // Advection
            fLow(fidx, loffset+lidx, delta2*Ap[lidx_f]*Vd*fac);

            fUpp(fidx, uoffset+uidx, (1-delta2)*delta1*Ap[lidx_f]*Vd*fac);
            if (delta1 != 1)
                fUpp(fidx, uoffset+uidx_m, (1-delta2)*(1-delta1)*Ap[lidx_f]*Vd*fac);

            //////////////////////
            // P-P diffusion (TODO)
            real_t dp = up[0]-lp[lnp-1];
            fLow(fidx, loffset+lidx, Dpp[lidx_f]*Vd*fac/dp);

            fUpp(fidx, uoffset+uidx, -delta1*Dpp[lidx_f]*Vd*fac/dp);
            if (delta1 != 1)
                fUpp(fidx, uoffset+uidx_m, -(1-delta1)*Dpp[lidx_f]*Vd*fac/dp);
        }

        loffset += lnp*lnxi;
        uoffset += unp*unxi;
    }
}
