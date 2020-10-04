/**
 * This boundary condition ensures that the trapped region is handled correctly
 * in kinetic simulations using toroidal geometry. The B.C. resets the equation
 * at f(xi0) corresponding -xi_T <= xi0 < 0 so that this region is exactly
 * equation to f(xi0) at 0 < -xi0 <= xi_T.
 */

#include <functional>
#include "FVM/Equation/BoundaryConditions/PXiInternalTrapping.hpp"


using namespace DREAM::FVM;
using namespace DREAM::FVM::BC;
using namespace std;


/**
 * Constructor.
 */
PXiInternalTrapping::PXiInternalTrapping(Grid *g)
    : BoundaryCondition(g) {
    
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

    this->trappedNegXi_indices = nullptr;
    this->nTrappedNegXi_indices = nullptr;
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
 *
 * This index array is passed to the PETSc function 'MatZeroRows()' when
 * resetting the equation.
 */
void PXiInternalTrapping::LocateTrappedRegion() {
    if (this->trappedNegXi_indices != nullptr)
        this->DeallocateTrappedIndices();

    len_t nr = this->grid->GetNr();
    this->trappedNegXi_indices = new PetscInt*[nr];
    this->trappedPosXi_indices = new PetscInt*[nr];
    this->nTrappedNegXi_indices = new PetscInt[nr];

    // Count number of xi0 cells inside [-xi_T, 0]
    for (len_t ir = 0; ir < nr; ir++) {
        auto mg = this->grid->GetMomentumGrid(ir);
        const real_t *xi0 = mg->GetP2();
        const real_t dxi0 = mg->GetDp2(0);
        const PetscInt nXi = mg->GetNp2();

        this->nTrappedNegXi_indices[ir] = 0;
        PetscInt minidx = nXi;

        for (PetscInt j = 0; j < nXi; j++) {
            if (this->grid->IsTrapped(ir, 0, j) && xi0[j] < 0) {
                this->nTrappedNegXi_indices[ir]++;

                if (j < minidx) minidx = j;
            }
        }

        // Populate index array (later to be passed to the PETSc Mat API)...
        this->trappedNegXi_indices[ir] = new PetscInt[this->nTrappedNegXi_indices[ir]];
        for (PetscInt i = 0; i < this->nTrappedNegXi_indices[ir]; i++)
            this->trappedNegXi_indices[ir][i] = minidx + i;

        // Find the closest mirrored xi0 (i.e. -xi0) on the grid...
        this->trappedPosXi_indices[ir] = new PetscInt[this->nTrappedNegXi_indices[ir]];
        for (PetscInt i = 0; i < this->nTrappedNegXi_indices[ir]; i++) {
            const len_t negXI0 = -xi0[this->trappedNegXi_indices[ir][i]];
            this->trappedPosXi_indices[ir][i] = nXi-1;

            for (PetscInt j = 0; j < nXi; j++) {
                // Allow for both monotonic increase and decrease in xi0...
                if ((dxi0 > 0 && xi0[j] > negXI0) || (dxi0 < 0 && xi0[j] < negXI0)) {
                    // Store the index right before (so we use j & j+1 for interpolation)...
                    if (j > 0)
                        j--;

                    this->trappedPosXi_indices[ir][i] = j;
                    break;
                }
            }
        }
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
void PXiInternalTrapping::AddToMatrixElements(Matrix*, real_t*) {
    // TODO
}

/**
 * Add elements to the function vector.
 */
void PXiInternalTrapping::AddToVectorElements(real_t*, const real_t*) {
    // TODO
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
    Matrix *mat, real_t*
) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    // Iterate over all rows with -xi_T <= xi0 < 0.
    for (len_t ir = 0; ir < nr; ir++) {
        PetscInt nidcs = this->nTrappedNegXi_indices[ir];
        PetscInt *idcs = this->trappedNegXi_indices[ir];

        // Clear the row
        mat->ZeroRows(nidcs, idcs);

        // NOTE: Must be 'INSERT_VALUES' since we called 'ZeroRows()' above
        // and didn't call 'PartialAssemble()' after...
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

        // Interpolation coefficient (xi0[J] is negative)
        // (TODO: This interpolation scheme most likely doesn't
        // preserve particle number, which we should require...)
        real_t delta =
            (pJ == nxi-1) ?
                1.0 :
                (-xi0[J] - xi0[pJ]) / (xi0[pJ+1] - xi0[pJ]);

        for (len_t i = 0; i < np; i++) {
            f(offset + J*np + i, offset + pJ*np + i, delta);

            if (pJ+1 < nxi-1)
                f(offset + J*np + i, offset + (pJ+1)*np + i, (1-delta));
        }
    }

    return np*nxi;
}

