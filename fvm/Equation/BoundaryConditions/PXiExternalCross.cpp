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

#include "FVM/Equation/BoundaryConditions/PXiExternalCross.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
PXiExternalCross::PXiExternalCross(
    DREAM::FVM::Grid *grid, DREAM::FVM::Grid *lowerGrid, DREAM::FVM::Grid *upperGrid,
    const DREAM::FVM::Operator*, enum condition_type ctype
) : BoundaryCondition(grid1), lowerGrid(lowerGrid), upperGrid(upperGrid),
    equation(eqn), type(ctype) {

    AllocateMemory();
}

/**
 * Destructor.
 */
PXiExternalCross::~PXiExternalCross() {
    DeallocateMemory();
}


bool PXiExternalCross::GridRebuild() {
    this->BoundaryCondition::GridRebuilt();

    this->DeallocateMemory();
    this->AllocateMemory();
    
    return true;
}

/**
 * Allocate memory for the interpolated advection/diffusion
 * coefficients.
 */
void PXiExternalCross::AllocateMemory() {
    // XXX here we assume that momentum grids are the same
    // at all radii
    const len_t nr  = this->upperGrid->GetNr();
    const len_t nxi = this->upperGrid->GetMomentumGrid(0)->GetNp2();

    this->upperAp  = new real_t[nr*nxi];
    this->upperDpp = new real_t[nr*nxi];
    this->upperDpx = new real_t[nr*nxi];
}

/**
 * Deallocate memory for interpolated coefficients.
 */
void PXiExternalCross::DeallocateMemory() {
    delete [] this->upperDpx;
    delete [] this->upperDpp;
    delete [] this->upperAp;
}

/**
 * Average and interpolate the advection/diffusion coefficients
 * on the lower grid onto the upper grid.
 */
bool PXiExternalCross::Rebuild(const real_t, UnknownQuantityHandler*) {
    const len_t nr = this->lowerGrid->GetNr();

    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *lmg = this->lowerGrid->GetMomentumGrid(ir);
        const MomentumGrid *umg = this->upperGrid->GetMomentumGrid(ir);

        const real_t *lxi  = lmg->GetP2();
        const real_t *uxi  = umg->GetP2();
        const real_t *lDxi = lmg->GetDp2();
        const real_t *uDxi = umg->GetDp2();

        const len_t lNp  = lmg->GetNp1();
        const len_t uNp  = umg->GetNp1();
        const len_t lNxi = lmg->GetNp2();
        const len_t uNxi = umg->GetNp2();

        const real_t *lVp_fp = this->lowerGrid->GetVp_f1(ir);
        const real_t *uVp_fp = this->upperGrid->GetVp_f1(ir);

        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        for (len_t j = 0, J = 0; j < uNxi; j++) {
            // Average coefficients...
            do {
                real_t d =
                    min(lxi_f[J+1], uxi_f[j+1]) - max(lxi_f[J], uxi_f[j]);
                real_t V =
                    lVp_fp[J*lNp + lNp] / uVp_fp[j*uNp + uNp];

                this->upperAp[ir*uNxi + j]  = V * (dxi/uDxi[j]) * Ap [J*lNp + lNp];
                this->upperDpp[ir*uNxi + j] = V * (dxi/uDxi[j]) * Dpp[J*lNp + lNp];
                this->upperDpx[ir*uNxi + j] = V * (dxi/uDxi[j]) * Dpx[J*lNp + lNp];

                if (uxi_f[j+1] > lxi_f[J+1] && J < lNxi-1)
                    J++;
                else
                    break;
            } while (true);
        }
    }

    return true;
}

/**
 * Add flux to jacobian block.
 */
void PXiExternalCross::AddToJacobianBlock(
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
void PXiExternalCross::AddToMatrixElements(
    Matrix *mat, real_t*
) {
    // Reset matrix offsets...
    const len_t rowOffset = mat->GetRowOffset();
    const len_t colOffset = mat->GetColOffset();
    mat->SetOffset(rowOffset, 0);

    /*#define f(I,J,V) mat->SetElement((I), (J), (V))
    #   include "PXiExternalCross.set.cpp"
    #undef f*/
    __SetElements([&mat](const len_t I, const len_t J, const real_t V) {
        mat->SetElement(I, J, V);
    });

    // Restore matrix offset
    mat->SetOffset(rowOffset, colOffset);
}

/**
 * Add flux to function vector.
 */
void PXiExternalCross::AddToVectorElements(
    real_t *vec, const real_t *f
) {
    /*#define f(I,J,V) vec[(I)] += (V)*f2[(J)]
    #   include "PXiExternalCross.set.cpp"
    #undef f*/
    __SetElements([&vec, &f](const len_t I, const len_t J, const real_t V) {
        vec[I] += V*f[J];
    });
}

/**
 * The flux across a boundary is given by the two components
 *
 *   Phi = F * ((1-d)*f_{Np} + d*f_{Np+1}),
 *
 * where 'f_{Np}' is the last point of the distribution on the
 * "lower" grid, and 'f_{Np+1}' is the first point of the distribution
 * on the "upper" grid. Depending on whether the ditribution function
 * is to be evaluated on the same grid as Phi or not, we either just
 * fetch its value, or calculate its average across the boundary.
 * This routine can evaluate both f_{Np} and f_{Np+1} (but only
 * evaluates one of them at a time).
 */
len_t PXiExternalCross::__SetFluxComponent(
    std::function<void(const len_t, const len_t, const real_t)>& f,
    const len_t firstIndex, const real_t coeff,
    const real_t xi_f1, const real_t xi_f2, const real_t dxi,
    const len_t lnxi, const real_t *lxi_f, const real_t *ldxi_f,
    const len_t jStart
) {
    // TODO
}

/**
 * Internal routine for setting matrix/vector elements.
 */
void PXiExternalCross::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> f
) {
    const len_t nr = this->grid->GetNr();
    len_t loffset = 0, uoffset = 0;

    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *lmg = this->lowerGrid->GetMomentumGrid(ir);
        const MomentumGrid *umg = this->upperGrid->GetMomentumGrid(ir);

        const len_t
            lnp = lmg->GetNp1(), lnxi = lmg->GetNp2(),
            unp = umg->GetNp2(), unxi = umg->GetNp2();

        const real_t
            *lxi_f = lmg->GetP2(),
            *uxi_f = umg->GetP2(),
            *ldp   = lmg->GetDp1(),
            *udp   = umg->GetDp1(),
            *ldxi  = lmg->GetDp2(),
            *udxi  = umg->GetDp2();

        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        const real_t
            *lVp   = this->lowerGrid->GetVp(ir),
            *uVp   = this->upperGrid->GetVp(ir);

        const real_t *delta1 = equation->GetInterpolationCoeff1(ir);

        if (this->type == TYPE_LOWER) {
            for (len_t j = 0, J = 0; j < lnxi; j++) {
                len_t idx  = j*lnp + lnp-1;
                len_t fidx = loffset + idx;

                real_t S = Ap[idx] / (lVp[idx] * ldp[lnp-1]);

                this->__SetFluxComponent(
                    f, fidx, -S*(1-delta1[idx]),
                    xi_f[j], xi_f[j+1], ldxi[j],
                    lnxi, lxi_f, ldxi_f, j
                );

                J = this->__SetFluxComponent(
                    f, fidx, S*(1-delta1[idx]),
                    xi_f[j], xi_f[j+1], ldxi,
                    unxi, uxi_f, J
                );
            }
        } else if (this->type == TYPE_UPPER) {
        } else if (this->type == TYPE_DENSITY) {
        } else
            throw EquationException("PXiExternalCross: Unrecognized boundary type: %d\n", this->type);

        /*for (len_t jl = 0, ju = 0; jl < lnxi; jl++) {
            const len_t idx1 = offset1 + j1*np1 +
                (this->type==TYPE_LOWER ? 0.0 : (np1-1));
            const len_t idx2 = offset2 + j2*np2 +
                (this->type==TYPE_LOWER ? (np2-1) : 0.0);

            // Index for coefficients
            const len_t cidx = (this->type==TYPE_LOWER ? idx2 : idx1);
            // Lower side of boundary
            const len_t lidx = (this->type==TYPE_LOWER ? idx2 : idx1);
            // Upper side of boundary
            const len_t uidx = (this->type==TYPE_LOWER ? idx1 : idx2);

            real_t S = 1.0 / (Vp[j1*np1+np1-1]*dp1[np1-1]*dxi1[j1]);

            // When we're considering the upper grid, we subtract
            // the flux from the lower grid.
            if (this->type == TYPE_UPPER)
                S = -S;

            // Patch together the flux across the boundary
            if (j2 == nxi2)
                j2--;

            do {
                real_t dxiBar =
                    min(xi1_f[j1+1], xi2_f[j2+1]) - max(xi1_f[j1], xi2_f[j1]);

                // Advection
                f(idx1, lidx, Ap[cidx]*(1-delta1[ir][cidx]) * S);
                f(idx1, uidx, Ap[cidx]*delta1[ir][cidx] * S);

                if (xi2_f[j2] < xi1_f[j1])
                    j2++;
            } while (j2 < nxi2 && xi2_f[j2] < xi1_f[j1]);
        }*/

        loffset += lnp*lnxi;
        uoffset += unp*unxi;
    }
}
