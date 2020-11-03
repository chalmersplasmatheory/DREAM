/**
 * This boundary condition ensures that the trapped region is handled correctly
 * in kinetic simulations using toroidal geometry. The B.C. resets the equation
 * at f(xi0) corresponding -xi_T <= xi0 < 0 so that this region is exactly
 * governed by the equation f(-xi0) = f(xi0) at 0 < -xi0 <= xi_T.
 * NOTE: Only the diagonal is reset, and the implementations of individual
 *       other equation terms must make sure that other elements are not set
 */

#include <functional>
#include "FVM/Equation/BoundaryConditions/PXiInternalTrapping.hpp"


using namespace DREAM::FVM;
using namespace DREAM::FVM::BC;
using namespace std;


/**
 * Constructor.
 */
PXiInternalTrapping::PXiInternalTrapping(
    Grid *g, Operator *oprtr
) : BoundaryCondition(g), fluxOperator(oprtr) {
    
    this->LocateTrappedRegion();
}

/**
 * Destructor.
 */
PXiInternalTrapping::~PXiInternalTrapping() {
    this->DeallocateTrappedIndices();
}

/**
 * Deallocate memory for the 'trappedNegXi_indices' array.
 */
void PXiInternalTrapping::DeallocateTrappedIndices() {
    len_t nr = this->grid->GetNr();
    for (len_t ir = 0; ir < nr; ir++)
        delete [] trappedNegXi_indices[ir];

    delete [] trappedNegXi_indices;
    delete [] nTrappedNegXi_indices;

    delete [] rowsToReset;

    this->trappedNegXi_indices = nullptr;
    this->nTrappedNegXi_indices = nullptr;
    this->rowsToReset = nullptr;
}

/**
 * Method which is called after the grid has been rebuilt.
 */
bool PXiInternalTrapping::GridRebuilt() {
    this->LocateTrappedRegion();
    return true;
}

/**
 * Construct the array 'trappedNegXi_indices', which contains indices to the
 * points of xi0 which are located inside the trapped region, and which have
 * xi0 < 0 (i.e. those points which should be copied from the points with
 * xi0 > 0).
 */
void PXiInternalTrapping::LocateTrappedRegion() {
    if (this->trappedNegXi_indices != nullptr)
        this->DeallocateTrappedIndices();

    len_t nr = this->grid->GetNr();
    this->trappedNegXi_indices = new PetscInt*[nr];
    this->trappedPosXi_indices = new PetscInt*[nr];
    this->nTrappedNegXi_indices = new PetscInt[nr];
    this->trappedNegXiRadial_indices = new PetscInt*[nr];
    this->trappedPosXiRadial_indices = new PetscInt*[nr];
    this->nTrappedNegXiRadial_indices = new PetscInt[nr];

    this->nRowsToReset = 0;
    // Count number of xi0 cells inside [-xi_T, 0]
    for (len_t ir = 0; ir < nr; ir++) {
        auto mg = this->grid->GetMomentumGrid(ir);
        const real_t *xi0   = mg->GetP2();
        const real_t *xi0_f = mg->GetP2_f();
        const real_t dxi0   = mg->GetDp2(0);
        const PetscInt nXi  = mg->GetNp2();
        const PetscInt nP   = mg->GetNp1();

        this->nTrappedNegXi_indices[ir] = 0;
        PetscInt minidx = nXi;
        PetscInt minidxRadial = nXi;

        for (PetscInt j = 0; j < nXi; j++){
            if (this->grid->IsNegativePitchTrappedIgnorableCell(ir,j)) {
                this->nTrappedNegXi_indices[ir]++;
                if (j < minidx)
                    minidx = j;
            }
            if (this->grid->IsNegativePitchTrappedIgnorableRadialFluxCell(ir,j)) {
                this->nTrappedNegXiRadial_indices[ir]++;
                if (j < minidxRadial) 
                    minidxRadial = j;
            }
        }

        // XI INDICES

        // Populate index array (later to be passed to the PETSc Mat API)...
        this->trappedNegXi_indices[ir] = new PetscInt[this->nTrappedNegXi_indices[ir]];
        for (PetscInt i = 0; i < this->nTrappedNegXi_indices[ir]; i++)
            this->trappedNegXi_indices[ir][i] = minidx + i;

        // Find the closest mirrored xi0 (i.e. -xi0) on the grid...
        this->trappedPosXi_indices[ir] = new PetscInt[this->nTrappedNegXi_indices[ir]];
        for (PetscInt i = 0; i < this->nTrappedNegXi_indices[ir]; i++) {
            const real_t negXI0 = -xi0[this->trappedNegXi_indices[ir][i]];
            this->trappedPosXi_indices[ir][i] = nXi-1;

            for (PetscInt j = 0; j < nXi; j++) {
                // Allow for both monotonic increase and decrease in xi0...
                if ((dxi0 > 0 && xi0_f[j+1] > negXI0) || (dxi0 < 0 && xi0_f[j+1] < negXI0)) {
                    this->trappedPosXi_indices[ir][i] = j;
                    break;
                }
            }
        }
        this->nRowsToReset += nP*nTrappedNegXi_indices[ir];

        // R INDICES
        // note: there is a lot of overlap here with the xi indices (sometimes the below is identical)
        //       but repeating it all here for now, for clarity

        // Populate index array (later to be passed to the PETSc Mat API)...
        this->trappedNegXiRadial_indices[ir] = new PetscInt[this->nTrappedNegXiRadial_indices[ir]];
        for (PetscInt i = 0; i < this->nTrappedNegXiRadial_indices[ir]; i++)
            this->trappedNegXiRadial_indices[ir][i] = minidxRadial + i;

        // Find the closest mirrored xi0 (i.e. -xi0) on the grid...
        this->trappedPosXiRadial_indices[ir] = new PetscInt[this->nTrappedNegXiRadial_indices[ir]];
        for (PetscInt i = 0; i < this->nTrappedNegXiRadial_indices[ir]; i++) {
            const real_t negXI0 = -xi0[this->trappedNegXiRadial_indices[ir][i]];
            this->trappedPosXiRadial_indices[ir][i] = nXi-1;

            for (PetscInt j = 0; j < nXi; j++) {
                // Allow for both monotonic increase and decrease in xi0...
                if ((dxi0 > 0 && xi0_f[j+1] > negXI0) || (dxi0 < 0 && xi0_f[j+1] < negXI0)) {
                    this->trappedPosXiRadial_indices[ir][i] = j;
                    break;
                }
            }
        }


    }

    // create large index vector for rows to reset 
    this->rowsToReset = new PetscInt[nRowsToReset];
    len_t offset=0;
    for (len_t ir = 0, it = 0; ir < nr; ir++) {
        len_t nP = grid->GetNp1(ir);
        len_t nXi = grid->GetNp2(ir);
        for(PetscInt j=0; j<nTrappedNegXi_indices[ir]; j++){
            len_t indXi = offset + nP*trappedNegXi_indices[ir][j];
            for (len_t i = 0; i < nP; i++, it++) 
                rowsToReset[it] = indXi + i;
        }
        offset += nP*nXi;
    }
}

/**
 * Rebuild this boundary condition (not used).
 */
bool PXiInternalTrapping::Rebuild(const real_t, UnknownQuantityHandler*) { return false; }

/**
 * Add elements to the Jacobian.
 */
void PXiInternalTrapping::AddToJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t*
) {
    if (uqtyId == derivId)
        this->AddToMatrixElements(jac, nullptr);
}

/**
 * Add elements to the linear operator matrix.
 */
void PXiInternalTrapping::AddToMatrixElements(Matrix *mat, real_t*) {
    this->_addElements([&mat](const len_t I, const len_t J, const real_t V) {
        mat->SetElement(I, J, V);
    });
}

/**
 * Add elements to the function vector.
 */
void PXiInternalTrapping::AddToVectorElements(real_t *vec, const real_t *f) {
    this->_addElements([&vec,&f](const len_t I, const len_t J, const real_t V) {
        vec[I] += V*f[J];
    });
}

/**
 * Internal routine for adding fluxes to matrix or vector. The purpose
 * of this routine is to add the flux from f(-xi_T^-}) to f(+xi_T^+),
 * i.e. from just _outside_ the trapped region at negative xi0, to just
 * _inside_ the trapped region at positive xi0.
 */
void PXiInternalTrapping::_addElements(
    function<void(const len_t, const len_t, const real_t)> f
) {
    const len_t nr = this->grid->GetNr();
    const enum AdvectionInterpolationCoefficient::adv_interp_mode interp_mode =
        AdvectionInterpolationCoefficient::AD_INTERP_MODE_FULL;

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        auto mg = this->grid->GetMomentumGrid(ir);
        const len_t np  = mg->GetNp1();
        const len_t nxi = mg->GetNp2();

        const real_t
            *Vp     = this->grid->GetVp(ir),
            *Vp_f2  = this->grid->GetVp_f2(ir),
            *dp_f   = mg->GetDp1(),
            *dxi0   = mg->GetDp2(),
            *dxi0_f = mg->GetDp2_f(),
            *dr     = this->grid->GetRadialGrid()->GetDr(),
            *dr_f   = this->grid->GetRadialGrid()->GetDr_f();

        const real_t *Ax  = fluxOperator->GetAdvectionCoeff2(ir);
        const real_t *Dxx = fluxOperator->GetDiffusionCoeff22(ir);
        const real_t *Dxp = fluxOperator->GetDiffusionCoeff21(ir);

        // indices indicating in which cells to mirror pitch fluxes
        len_t jm = trappedNegXi_indices[ir][0];
        len_t jp = trappedPosXi_indices[ir][0];

        // Iterate over all p and set the fluxes...
        for (len_t i = 0; i < np; i++) {
            const len_t idxm = jm*np + i;
            const len_t idxp = jp*np + i;

            // XI ADVECTION
            real_t S_i = Ax[idxm] * Vp_f2[idxm] / (Vp[idxp]*dxi0[jp]);
            AdvectionInterpolationCoefficient *delta2 = fluxOperator->GetInterpolationCoeff2();
            const real_t *delta = delta2->GetCoefficient(ir, i, jm, interp_mode);
            for (len_t n, k = delta2->GetKmin(jm, &n); k <= delta2->GetKmax(jm, nxi); k++, n++)
                f(offset+idxp, offset+k*np+i, -S_i * delta[n]);

            // XI-XI DIFFUSION
            if (jm > 0) {
                S_i = Dxx[idxm] * Vp_f2[idxm] / (Vp[idxp]*dxi0[jp]*dxi0_f[jm-1]);

                f(offset+idxp, offset+jm*np+i,     +S_i);
                f(offset+idxp, offset+(jm-1)*np+i, -S_i);
            }

            // XI-P DIFFUSION
            if (jm > 0 && (i > 0 && i < np-1)) {
                S_i = Dxp[idxm] * Vp_f2[idxm] / (Vp[idxp]*dxi0[jp]*(dp_f[i]+dp_f[i-1]));

                f(offset+idxp, offset+(jm-1)*np+i+1, +S_i);
                f(offset+idxp, offset+jm*np+i+1,     +S_i);
                f(offset+idxp, offset+(jm-1)*np+i-1, -S_i);
                f(offset+idxp, offset+jm*np+i-1,     -S_i);
            }
        }

        const real_t *Ar  = fluxOperator->GetAdvectionCoeffR(ir);
        const real_t *Drr = fluxOperator->GetDiffusionCoeffRR(ir);

        // unlike the xi fluxes, here there may be multiple radial fluxes to move for 
        // each radius, so we sum over all pitch cells that may be mirrored
        for(len_t j=0; j<(len_t) nTrappedNegXiRadial_indices[ir]; j++){
            // pitch indices
            len_t km = trappedNegXiRadial_indices[ir][j];
            len_t kp = trappedPosXiRadial_indices[ir][j];
            
            // Iterate over all p and set the fluxes...
            for (len_t i = 0; i < np; i++) {
                const len_t idxm = km*np + i;
                const len_t idxp = kp*np + i;

                // R ADVECTION
                real_t S_i = Ar[idxm] * Vp_f2[idxm] / (Vp[idxp]*dr[ir]);
                if(S_i){ // often we will not have radial fluxes, and can skip these calculations
                    AdvectionInterpolationCoefficient *deltar = fluxOperator->GetInterpolationCoeffR();
                    const real_t *delta = deltar->GetCoefficient(ir, i, km, interp_mode);
                    for (len_t n, k = deltar->GetKmin(ir, &n); k <= deltar->GetKmax(ir, nr); k++, n++)
                        f(offset+idxp, offset+(k-ir)*np*nxi + idxm, -S_i * delta[n]);
                }

                // R-R DIFFUSION
                if (ir > 0) {
                    S_i = Drr[idxm] * Vp_f2[idxm] / (Vp[idxp]*dr[ir]*dr_f[ir-1]);
                    if(S_i){
                        f(offset+idxp, offset + idxm,        +S_i);
                        f(offset+idxp, offset-np*nxi + idxm, -S_i);
                    }
                }
            }
        }

        offset += np*nxi;
    }
}

/**
 * Reset the rows of the jacobian corresponding to -xi_T <= xi0 < 0.
 */
void PXiInternalTrapping::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t*
) {
    if (uqtyId == derivId)
        this->SetMatrixElements(jac, nullptr);
}

/**
 * Reset the rows of the matrix corresponding to -xi_T <= xi0 < 0
 * and set f(xi0) = f(-xi0) in this region.
 */
void PXiInternalTrapping::SetMatrixElements(
    Matrix *mat, real_t *rhs
) {
    // override diagonal and reset rhs in trapped region
    mat->SetDiagonalConstant(nRowsToReset, rowsToReset, -1.0);
    if(rhs != nullptr)
        for(len_t it=0; it<nRowsToReset; it++)
            rhs[rowsToReset[it]] = 0;
    
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;
    // Iterate over all rows with -xi_T <= xi0 < 0.
    for (len_t ir = 0; ir < nr; ir++) {
        // NOTE: Must be 'INSERT_VALUES' since we called 'SetDiagonalConstant(...)' 
        // above and didn't call 'PartialAssemble()' after...
        offset += this->_setElements(
            ir, offset,
            [&mat](const len_t I, const len_t J, const real_t v)
            { mat->SetElement(I, J, v, INSERT_VALUES); }
        );
    }
}

/**
 * Reset the elements of the vector corresponding to -xi_T <= xi0 < 0
 * and set f(xi0) = f(-xi0) in this region.
 */
void PXiInternalTrapping::SetVectorElements(real_t *vec, const real_t *f) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    // Iterate over all rows with -xi_T <= xi0 < 0.
    for (len_t ir = 0; ir < nr; ir++) {
        offset += this->_setElements(
            ir, offset,
            [&vec,&f](const len_t I, const len_t J, const real_t v) { vec[I] = v*f[J]; }
        );
    }
}

/**
 * Internal routine used for setting f(xi0) = f(-xi0) in the region
 * xi_T <= xi0 < 0.
 *
 * offset: Offset of the current radial coordinate.
 * f:      Function which sets the element.
 */
len_t PXiInternalTrapping::_setElements(
    const len_t ir, const len_t offset,
    function<void(const len_t, const len_t, const real_t)> f
) {
    PetscInt nidcs = this->nTrappedNegXi_indices[ir];
    PetscInt *idcs = this->trappedNegXi_indices[ir];
    // Set the row to f(xi0) = f(-xi0)
    auto mg = this->grid->GetMomentumGrid(ir);
    const len_t np  = mg->GetNp1();
    const len_t nxi = mg->GetNp2();
    const real_t *xi0 = mg->GetP2();
    
    for (PetscInt j = 0; j < nidcs; j++) {
        const len_t J  = idcs[j];
        const len_t pJ = this->trappedPosXi_indices[ir][j];

        // interpolate in the direction of the closest grid point
        int_t interpolationDirection;
        if( xi0[pJ] > - xi0[J])
            interpolationDirection = -1;
        else 
            interpolationDirection = 1;

        // Interpolation coefficient (xi0[J] is negative)
        real_t delta =
            (pJ == nxi-1) ?
                1.0 :
                (-xi0[J] - xi0[pJ]) / (xi0[pJ+interpolationDirection] - xi0[pJ]);
        // set nearly vanishing elements to identically zero to reduce nnz
        if(abs(1-delta)<100*realeps)
            delta=1.0;
        else if (abs(delta)<100*realeps)
            delta=0.0;


        for (len_t i = 0; i < np; i++) {
            f(offset + J*np + i, offset + pJ*np + i, 1-delta);

            if (pJ+1 < nxi-1)
                f(offset + J*np + i, offset + (pJ+interpolationDirection)*np + i, delta);
        }
    }

    return np*nxi;
}
