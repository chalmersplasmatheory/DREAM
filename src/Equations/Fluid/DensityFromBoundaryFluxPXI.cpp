/**
 * This term represents the flux of particles from a p/xi grid across the grid
 * p boundary into the fluid grid. The flux of particles into the fluid grid
 * cell 'r' is given by the pitch-integrated flux across the p=pmax boundary:
 *
 *                      / +1
 *   S(r) = 2*pi*pmax^2 |    Phi^{(p)}(r, pmax, xi) dxi
 *                      / -1
 *
 * NOTE: This term should be added to the equation for the density (e.g. 'nRE')
 * and applied to the distribution function (i.e. the matrix column representing 'fRE').
 */

#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/NotImplementedException.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
DensityFromBoundaryFluxPXI::DensityFromBoundaryFluxPXI(
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid,
    const FVM::Equation *eqn
) : FVM::EquationTerm(densityGrid), distributionGrid(distributionGrid), equation(eqn) { }


/**
 * Destructor.
 */
DensityFromBoundaryFluxPXI::~DensityFromBoundaryFluxPXI() {}


/**
 * Get the number of non-zero elements per matrix row in the
 * linearized operator matrix.
 */
len_t DensityFromBoundaryFluxPXI::GetNumberOfNonZerosPerRow() const {
    // XXX Here we assume that all momentum grids are the same
    const len_t nXi = this->distributionGrid->GetMomentumGrid(0)->GetNp2();

    return (nXi*3);
}

/**
 * Get the number of non-zero elements per matrix row in the
 * jacobian matrix.
 *
 * TODO calculate correctly by taking 'equation' into account.
 */
len_t DensityFromBoundaryFluxPXI::GetNumberOfNonZerosPerRow_jac() const {
    // XXX Here we assume that all momentum grids are the same
    const len_t nXi = this->distributionGrid->GetMomentumGrid(0)->GetNp2();

    return (nXi*3);
}

/**
 * Set a block of the jacobian.
 *
 * derivId: Quantity with respect to which the differentation should be made.
 * qtyId:   Quantity to differentiate.
 * jac:     Jacobian matrix.
 */
void DensityFromBoundaryFluxPXI::SetJacobianBlock(
    const len_t /*derivId*/, const len_t /*qtyId*/, FVM::Matrix * /*jac*/
) {
    throw NotImplementedException("Cannot set jacobian for 'DensityFromBoundaryFluxPXI' term yet.");
}

/**
 * Set elements in the linearized operator matrix.
 *
 * mat: Matrix to set elements of.
 * rhs: Right-hand-side vector to set elements of.
 */
void DensityFromBoundaryFluxPXI::SetMatrixElements(
    FVM::Matrix *mat, real_t*
) {
    const len_t nr = this->grid->GetNr();
    const real_t *VpVol = this->grid->GetVpVol();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = this->grid->GetMomentumGrid(ir);
        const len_t
            nxi = mg->GetNp2(),
            np  = mg->GetNp1();

        const real_t
            *dp    = mg->GetDp1(),
            *dxi   = mg->GetDp2(),
            *dp_f  = mg->GetDp1_f(),
            *Vp_fp = this->grid->GetVp_f1(ir);

        const real_t *Ap = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        // Evaluate xi-integral
        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = offset + j*np + (np-1);

            // Contribution from advection (with delta = 0)
            mat->SetElement(ir, idx, Ap[j*(np+1) + np]);

            // Constribution from diffusion
            real_t dd = dp[np-1] / dp[np-2];
            real_t dVol = Vp_fp[j*(np+1) + np] * dxi[j] / VpVol[ir];

            // Dpp
            mat->SetElement(ir, idx, (1+dd) * Dpp[j*(np+1) + np-1] / dp_f[np-2] * dVol);
            mat->SetElement(ir, idx-1, -(1+dd) * Dpp[j*(np+1) + np-1] / dp_f[np-2] * dVol);

            mat->SetElement(ir, idx-1, -dd * Dpp[j*(np+1) + np-2] / dp_f[np-3] * dVol);
            mat->SetElement(ir, idx-2,  dd * Dpp[j*(np+1) + np-2] / dp_f[np-3] * dVol);

            // TODO Dpx
        }

        offset += np*nxi;
    }
}

/**
 * Set elements in the given function vector.
 *
 * vec: Function vector to set elements of.
 * f:   Distribution function.
 */
void DensityFromBoundaryFluxPXI::SetVectorElements(
    real_t *vec, const real_t *f
) {
    const len_t nr = this->grid->GetNr();
    const real_t *VpVol = this->grid->GetVpVol();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = this->grid->GetMomentumGrid(ir);
        const len_t
            nxi = mg->GetNp2(),
            np  = mg->GetNp1();

        const real_t
            *dp    = mg->GetDp1(),
            *dxi   = mg->GetDp2(),
            *dp_f  = mg->GetDp1_f(),
            *Vp_fp = this->grid->GetVp_f1(ir);

        len_t Sr = 0;
        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        // Evaluate xi-integral
        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = offset + j*np + (np-1);

            // Contribution from advection (with delta = 0)
            real_t PhiP_a = Ap[j*(np+1) + np] * f[idx];

            // Constribution from diffusion
            // Dpp
            real_t Phi12  = Dpp[j*(np+1) + np-1] * (f[idx] - f[idx-1]) / dp_f[np-2];
            real_t Phi32  = Dpp[j*(np+1) + np-2] * (f[idx-1] - f[idx-2]) / dp_f[np-3];
            real_t PhiP_d = Phi12 + dp[np-1] / dp[np-2] * (Phi12 - Phi32);

            // TODO Dpx

            Sr += (PhiP_a + PhiP_d) * Vp_fp[j*(np+1) + np] * dxi[j];
        }

        vec[ir] += Sr / VpVol[ir];

        offset += np*nxi;
    }
}

