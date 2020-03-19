/**
 * Implementation of the internal boundary condition on the xi = +-1 boundaries.
 * Due to symmetry, the flux across this boundary must vanish, and hence this
 * boundary condition is trivial and doesn't need to set any matrix elements.
 */

#include "FVM/Grid/PXiGrid/PInternalBoundaryCondition.hpp"

using namespace TQS::FVM::BC;


/**
 * Allocate memory for the enforced momentum-space fluxes.
 */
void PInternalBoundaryCondition::AllocateFluxes() {
    nxi = new len_t[nr];
    p2S = new real_t*[nr];

    for (len_t i = 0; i < nr; i++) {
        nxi[i] = grid->GetMomentumGrid(i)->GetNp1();
        p2S[i] = new real_t[nxi[i]];
    }
}

/**
 * Deallocate memory for the enforced momentum-space fluxes.
 */
void PInternalBoundaryCondition::DeallocateFluxes() {
    for (len_t i = 0; i < nr; i++)
        delete [] p2S[i];

    delete [] p2S;
    delete [] nxi;
}

/**
 * Function called when the underlying 'RadialGrid' is re-built.
 */
bool PInternalBoundaryCondition::GridRebuilt() {
    DeallocateFluxes();

    nr  = grid->GetNr();

    AllocateFluxes();
}

/**
 * Function called when coefficients for operators should be
 * re-built. Since this boundary condition always the same
 * (zero flux across xi = +-1), we don't need to do anything
 * special here.
 */
bool PInternalBoundaryCondition::Rebuild(const real_t) { return false; }

/**
 * Set the matrix elements corresponding to this boundary
 * condition. Since the condition is trivial (zero flux across xi = +-1)
 * we don't need to set any matrix element explicitly.
 *
 * mat: Matrix to set elements in.
 */
void PInternalBoundaryCondition::SetMatrixElements(Matrix*, real_t *rhs) {
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const len_t nxi = this->nxi[ir];
        const len_t np  = grid->GetMomentumGrid(ir)->GetNp1();

        for (len_t j = 0; j < nxi; j++) {
            // Modify RHS vector
            rhs[offset + j*np] += this->p2S[ir][j];
        }

        offset += np * nxi;
    }
}

