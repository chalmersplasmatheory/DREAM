/**
 * An external boundary condition which represents an open external boundary
 * which permits outflux of particles. This boundary condition is designed to
 * be used in conjunction with the 'DensityFromBoundaryFluxPXi' equation term,
 * which adds particles to a fluid quantity 'n' at the same rate as they are
 * lost across the boundary handled by this boundary condition.
 */

#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"


using namespace DREAM::FVM::BC;


/**
 * Constructor.
 */
PXiExternalLoss::PXiExternalLoss(Grid *g, const Equation *eqn)
    : BoundaryCondition(g), equation(eqn) { }

/**
 * Destructor.
 */
PXiExternalLoss::~PXiExternalLoss() { }


/**
 * Rebuild coefficients for this term.
 * (not used)
 */
bool PXiExternalLoss::Rebuild(const real_t, UnknownQuantityHandler*) { return false; }

/**
 * Add flux to jacobian block.
 */
void PXiExternalLoss::AddToJacobianBlock(
    const len_t /*derivId*/, const len_t /*qtyId*/, Matrix * /*jac*/
) {
    throw EquationException("PXiExternalLoss: AddToJacobianBlock(): Not implemented yet.");
}

/**
 * Add flux to linearized operator matrix.
 *
 * mat: Matrix to add boundary conditions to.
 * rhs: Right-hand-side vector (not used).
 */
void PXiExternalLoss::AddToMatrixElements(
    Matrix *mat, real_t*
) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;
    
    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *mg = this->grid->GetMomentumGrid(ir);
        const len_t
            np  = mg->GetNp1(),
            nxi = mg->GetNp2();

        const real_t
            *dp    = mg->GetDp1(),
            *dp_f  = mg->GetDp1_f(),
            *Vp    = this->grid->GetVp(ir),
            *Vp_fp = this->grid->GetVp_f1(ir);

        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = offset + j*np + (np-1);
            const real_t Vd = Vp_fp[j*(np+1) + np] / (Vp[idx-offset] * dp[np-1]);
            const real_t dd = dp[np-1] / dp[np-2];

            // Contribution from advection (with delta = 0)
            mat->SetElement(idx, idx, -Ap[j*(np+1) + np] * Vd);

            // Contribution from diffusion
            // Dpp
            mat->SetElement(idx, idx,   -(1+dd)*Dpp[j*(np+1) + np-1]/dp_f[np-2] * Vd);
            mat->SetElement(idx, idx-1,  (1+dd)*Dpp[j*(np+1) + np-1]/dp_f[np-2] * Vd);

            mat->SetElement(idx, idx-1,  dd*Dpp[j*(np+1) + np-2]/dp_f[np-3] * Vd);
            mat->SetElement(idx, idx-2, -dd*Dpp[j*(np+1) + np-2]/dp_f[np-3] * Vd);

            // Dpx
            if (j > 0 && j < nxi-1) {
                mat->SetElement(idx, idx+np,   -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);
                mat->SetElement(idx, idx+np-1, -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);
                mat->SetElement(idx, idx-np,   +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);
                mat->SetElement(idx, idx-np-1, +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);

                mat->SetElement(idx, idx+np-1, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
                mat->SetElement(idx, idx+np-2, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
                mat->SetElement(idx, idx-np-1, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
                mat->SetElement(idx, idx-np-2, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
            }
        }

        offset += np*nxi;
    }
}

/**
 * Add flux to function vector.
 *
 * vec: Function vector to add boundary conditions to.
 * f:   Current value of distribution function.
 */
void PXiExternalLoss::AddToVectorElements(
    real_t *vec, const real_t *f
) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;
    
    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *mg = this->grid->GetMomentumGrid(ir);
        const len_t
            np  = mg->GetNp1(),
            nxi = mg->GetNp2();

        const real_t
            *dp    = mg->GetDp1(),
            *dp_f  = mg->GetDp1_f(),
            *Vp    = this->grid->GetVp(ir),
            *Vp_fp = this->grid->GetVp_f1(ir);

        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        for (len_t j = 0; j < nxi; j++) {
            const len_t idx = offset + j*np + (np-1);

            // Contribution from advection (with delta = 0)
            real_t PhiP_a = Ap[j*(np+1) + np] * f[idx];

            // Contribution from diffusion
            // Dpp
            real_t Phi12  = Dpp[j*(np+1) + np-1] * (f[idx]- f[idx-1]) / dp_f[np-2];
            real_t Phi32  = Dpp[j*(np+1) + np-2] * (f[idx-1] - f[idx-2]) / dp_f[np-3];

            // Dpx
            if (j > 0 && j < nxi-1) {
                Phi12 += Dpx[j*(np+1) + np-1] * (f[idx+np+1]+f[idx+np] - f[idx-np+1]-f[idx-np]) / (dp_f[np-2]+dp_f[np-3]);
                Phi32 += Dpx[j*(np+1) + np-2] * (f[idx+np]+f[idx+np-1] - f[idx-np]-f[idx-np-1]) / (dp_f[np-3]+dp_f[np-4]);
            }

            real_t PhiP_d = Phi12 + dp[np-1] / dp[np-2] * (Phi12 - Phi32);

            real_t PhiP = PhiP_a + PhiP_d;

            vec[idx] -= PhiP * Vp_fp[j*(np+1) + np] / (Vp[idx-offset] * dp[np-1]);
        }

        offset += np*nxi;
    }
}

