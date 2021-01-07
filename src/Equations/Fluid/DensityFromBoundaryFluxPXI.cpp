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

#include <functional>
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/NotImplementedException.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
DensityFromBoundaryFluxPXI::DensityFromBoundaryFluxPXI(
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid,
    const FVM::Operator *eqn
) : FVM::EquationTerm(densityGrid), distributionGrid(distributionGrid), 
    equation(eqn) {}


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
    if (qtyId == derivId)
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
    this->__SetElements([&mat](const len_t I, const len_t J, const real_t V) {
        mat->SetElement(I, J, V);
    });
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
    this->__SetElements([&vec,&f](const len_t I, const len_t J, const real_t V) {
        vec[I] += V*f[J];
    });
}

void DensityFromBoundaryFluxPXI::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> f
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

        // Evaluate xi-integral
        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = offset + j*np + (np-1);

            //real_t dd = dp[np-1] / dp[np-2];
            real_t dd = 0;
            real_t dVol = Vp_fp[j*(np+1) + np-1] * dxi[j] / VpVol[ir];

            // Contribution from advection (with delta = 0)
            const real_t *delta1_0 = equation->GetInterpolationCoeff1(ir,np-1,j);
            for(len_t k=0; k<3; k++)
                f(ir, idx+k-2, -(1+dd)*delta1_0[k]*Ap[j*(np+1) + np-1] * dVol);

//            f(ir, idx, -(1+dd)*delta0[2]*Ap[j*(np+1)+np-1] * dVol);
//            f(ir, idx-1, -(1+dd)*delta0[1]*Ap[j*(np+1)+np-1] * dVol);

            const real_t *delta1_1 = equation->GetInterpolationCoeff1(ir,np-2,j);
            for(len_t k=0; k<4; k++)
                f(ir, idx+k-3, dd*delta1_1[k]*Ap[j*(np+1) + np-2] * dVol);
//            f(ir, idx-1, dd*delta1[2]*Ap[j*(np+1)+np-2] * dVol);
//            f(ir, idx-2, dd*delta1[1]*Ap[j*(np+1)+np-2] * dVol);

            // Constribution from diffusion
            // Dpp
            f(ir, idx,    (1+dd) * Dpp[j*(np+1) + np-1] / dp_f[np-2] * dVol);
            f(ir, idx-1, -(1+dd) * Dpp[j*(np+1) + np-1] / dp_f[np-2] * dVol);

            f(ir, idx-1, -dd * Dpp[j*(np+1) + np-2] / dp_f[np-3] * dVol);
            f(ir, idx-2,  dd * Dpp[j*(np+1) + np-2] / dp_f[np-3] * dVol);

            // Dpx
            if (j > 0 && j < nxi-1) {
                f(ir, idx+np,   +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);
                f(ir, idx+np-1, +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);
                f(ir, idx-np,   -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);
                f(ir, idx-np-1, -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * dVol);

                f(ir, idx+np-1, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
                f(ir, idx+np-2, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
                f(ir, idx-np-1, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
                f(ir, idx-np-2, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * dVol);
            }
        }

        offset += np*nxi;
    }
}

