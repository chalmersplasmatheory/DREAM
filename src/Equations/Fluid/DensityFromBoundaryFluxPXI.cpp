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
    const FVM::Equation *eqn, len_t fId, len_t momentId
) : FVM::EquationTerm(densityGrid), distributionGrid(distributionGrid), 
    equation(eqn), fId(fId), momentId(momentId) { }


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
 * x:       Value of unknown quantity.
 */
void DensityFromBoundaryFluxPXI::SetJacobianBlock(
    const len_t qtyId, const len_t derivId, FVM::Matrix * jac, const real_t* /*x*/
) {
    //throw NotImplementedException("Cannot set jacobian for 'DensityFromBoundaryFluxPXI' term yet.");
    if ((qtyId==momentId) && (derivId == fId))
        this->SetMatrixElements(jac, nullptr);

    // TODO Handle derivatives of coefficients
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
    const len_t nr = this->distributionGrid->GetNr();
    const real_t *VpVol = this->distributionGrid->GetVpVol();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = this->distributionGrid->GetMomentumGrid(ir);
        const len_t
            nxi = mg->GetNp2(),
            np  = mg->GetNp1();

        const real_t
            //*dp    = mg->GetDp1(),
            *dxi   = mg->GetDp2(),
            *dp_f  = mg->GetDp1_f(),
            *Vp_fp = this->distributionGrid->GetVp_f1(ir);

        const real_t *Ap = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        const real_t *delta1 = equation->GetInterpolationCoeff1(ir);

        // Evaluate xi-integral
        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = offset + j*np + (np-1);

            //real_t dd = dp[np-1] / dp[np-2];
            real_t dd = 0;
            real_t dVol = Vp_fp[j*(np+1) + np-1] * dxi[j] / VpVol[ir];

            // Contribution from advection (with delta = 0)
            mat->SetElement(ir, idx, -(1+dd)*delta1[j*np+(np-1)]*Ap[j*(np+1)+np-1] * dVol);
            mat->SetElement(ir, idx-1, -(1+dd)*(1-delta1[j*np+(np-1)])*Ap[j*(np+1)+np-1] * dVol);

            mat->SetElement(ir, idx-1, dd*delta1[j*np+(np-2)]*Ap[j*(np+1)+np-2] * dVol);
            mat->SetElement(ir, idx-2, dd*(1-delta1[j*np+(np-2)])*Ap[j*(np+1)+np-2] * dVol);

            // Constribution from diffusion
            // Dpp
            mat->SetElement(ir, idx,    (1+dd) * Dpp[j*(np+1) + np-1] / dp_f[np-2] * dVol);
            mat->SetElement(ir, idx-1, -(1+dd) * Dpp[j*(np+1) + np-1] / dp_f[np-2] * dVol);

            mat->SetElement(ir, idx-1, -dd * Dpp[j*(np+1) + np-2] / dp_f[np-3] * dVol);
            mat->SetElement(ir, idx-2,  dd * Dpp[j*(np+1) + np-2] / dp_f[np-3] * dVol);

            // Dpx
            if (j > 0 && j < nxi-1) {
                mat->SetElement(ir, idx+np,   +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);
                mat->SetElement(ir, idx+np-1, +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);
                mat->SetElement(ir, idx-np,   -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);
                mat->SetElement(ir, idx-np-1, -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);

                mat->SetElement(ir, idx+np-1, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
                mat->SetElement(ir, idx+np-2, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
                mat->SetElement(ir, idx-np-1, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
                mat->SetElement(ir, idx-np-2, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
            }
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
    const len_t nr = this->distributionGrid->GetNr();
    const real_t *VpVol = this->distributionGrid->GetVpVol();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = this->distributionGrid->GetMomentumGrid(ir);
        const len_t
            nxi = mg->GetNp2(),
            np  = mg->GetNp1();

        const real_t
            //*dp    = mg->GetDp1(),
            *dxi   = mg->GetDp2(),
            *dp_f  = mg->GetDp1_f(),
            *Vp_fp = this->distributionGrid->GetVp_f1(ir);

        real_t Sr = 0;
        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        const real_t *delta1 = equation->GetInterpolationCoeff1(ir);

        // Evaluate xi-integral
        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = offset + j*np + (np-1);

            // Contribution from advection (with delta = 0)
            //real_t PhiP_a = Ap[j*(np+1) + np] * f[idx];
            real_t PhiP_a12 = Ap[j*(np+1) + np-1] * (
                delta1[j*np+(np-1)]*f[idx] + (1-delta1[j*np+(np-1)])*f[idx-1]
            );
            real_t PhiP_a32 = Ap[j*(np+1) + np-2] * (
                delta1[j*np+(np-2)]*f[idx-1] + (1-delta1[j*np+(np-2)])*f[idx-2]
            );

            // Constribution from diffusion
            // Dpp
            real_t PhiP_d12  = Dpp[j*(np+1) + np-1] * (f[idx] - f[idx-1]) / dp_f[np-2];
            real_t PhiP_d32  = Dpp[j*(np+1) + np-2] * (f[idx-1] - f[idx-2]) / dp_f[np-3];

            // Dpx
            if (j > 0 && j < nxi-1) {
                PhiP_d12 += Dpx[j*(np+1) + np-1] * (f[idx+np]+f[idx+np-1] - f[idx-np]-f[idx-np-1]) / (dp_f[np-2]+dp_f[np-3]);
                PhiP_d32 += Dpx[j*(np+1) + np-2] * (f[idx+np-1]+f[idx+np-2] - f[idx-np-1]-f[idx-np-2]) / (dp_f[np-3]+dp_f[np-4]);
            }

            //const real_t dd = dp[np-1] / dp[np-2];
            const real_t dd = 0;

            real_t PhiP_a = PhiP_a12 + dd * (PhiP_a12 - PhiP_a32);
            real_t PhiP_d = PhiP_d12 + dd * (PhiP_d12 - PhiP_d32);

            real_t PhiP = -PhiP_a + PhiP_d;

            Sr += PhiP * Vp_fp[j*(np+1) + np-1] * dxi[j];
        }

        vec[ir] += Sr / VpVol[ir];

        offset += np*nxi;
    }
}

