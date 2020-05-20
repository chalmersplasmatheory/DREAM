/**
 * Implementation of the internal boundary condition on the xi = +-1 boundaries.
 * Due to symmetry, the flux across this boundary must vanish, and hence this
 * boundary condition is trivial and doesn't need to set any matrix elements.
 */

#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"

using namespace DREAM::FVM::BC;


PInternalBoundaryCondition::PInternalBoundaryCondition(Grid *g)
    : BoundaryCondition(g) {

    this->nr  = g->GetNr();
    
    this->AllocateFluxes();
}

/**
 * Allocate memory for the enforced momentum-space fluxes.
 */
void PInternalBoundaryCondition::AllocateFluxes() {
    nxi = new len_t[nr];
    VpS = new real_t*[nr];

    for (len_t i = 0; i < nr; i++) {
        nxi[i] = grid->GetMomentumGrid(i)->GetNp2();
        VpS[i] = new real_t[nxi[i]];
    }

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t nxi = this->nxi[ir];

        for (len_t j = 0; j < nxi; j++)
            VpS[ir][j] = 0;
    }
}

/**
 * Deallocate memory for the enforced momentum-space fluxes.
 */
void PInternalBoundaryCondition::DeallocateFluxes() {
    for (len_t i = 0; i < nr; i++)
        delete [] VpS[i];

    delete [] VpS;
    delete [] nxi;
}

/**
 * Function called when the underlying 'RadialGrid' is re-built.
 */
bool PInternalBoundaryCondition::GridRebuilt() {
    DeallocateFluxes();

    this->nr = grid->GetNr();

    AllocateFluxes();
    return true;
}

/**
 * Function called when coefficients for operators should be
 * re-built. Since this boundary condition always the same
 * (zero flux across xi = +-1), we don't need to do anything
 * special here.
 */
bool PInternalBoundaryCondition::Rebuild(const real_t, UnknownQuantityHandler*) { return false; }

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
            rhs[offset + j*np] += this->VpS[ir][j];
        }

        offset += np * nxi;
    }
}

