/**
 * An external boundary condition which represents an open external boundary
 * which permits outflux of particles. This boundary condition is designed to
 * be used in conjunction with the 'DensityFromBoundaryFluxPXi' equation term,
 * which adds particles to a fluid quantity 'n' at the same rate as they are
 * lost across the boundary handled by this boundary condition.
 */

#include <functional>
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"


using namespace DREAM::FVM::BC;


/**
 * Constructor.
 *
 * g:        Grid on which the unknown quantity that this quantity is applied
 *           to lives.
 * eqn:      Operator which causes the flux of particles across this boundary.
 * fId:      ID of the distribution function involved in this B.C.
 * targetId: ID of the unknown quantity which is the target of this B.C.
 *           I.e., if this B.C. acts as a source on fluid grid, then
 *           'targetId' should be the ID of the fluid quantity. Otherwise,
 *           it should be the same as 'fId'.
 * distGrid: If 'boundary' is 'BOUNDARY_FLUID', then 'g' is the fluid grid,
 *           and one has to give the distribution function grid here. If
 *           'nullptr', then 'g' is assumed to be the distribution grid.
 * bc:       Type of boundary condition to represent.
 * boundary: Which side of the boundary this condition applies to.
 */
PXiExternalLoss::PXiExternalLoss(
    Grid *g, const Operator *eqn, const len_t fId, const len_t targetId,
    Grid *distGrid, 
    enum boundary_type boundary, enum bc_type bc
) : BoundaryCondition(g), equation(eqn), fId(fId), targetId(targetId),
    boundaryCondition(bc), boundary(boundary) {
    
    if (distGrid == nullptr) {
        this->distributionGrid = g;
        
        if (boundary != BOUNDARY_KINETIC)
            throw OperatorException(
                "Invalid configuration of boundary condition 'PXiExternalLoss'. "
                "Boundary is of type 'fluid', but no kinetic grid was given."
            );
    } else
        this->distributionGrid = distGrid;
}

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
    const len_t qtyId, const len_t derivId, Matrix * jac, const real_t* /*x*/
) {
    if ((derivId == this->fId) && (this->fId == qtyId))
        this->AddToMatrixElements(jac, nullptr);

    // TODO handle derivatives of coefficients
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
    this->__SetElements([&mat](const len_t I, const len_t J, const real_t V) {
        mat->SetElement(I, J, V);
    });
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
    this->__SetElements([&vec,&f](const len_t I, const len_t J, const real_t V) {
        vec[I] += V*f[J];
    });
}

/**
 * PRIVATE
 * Internal routine used for setting matrix/vector elements.
 *
 * f(I,J,V): Function for setting matrix/vector elements. I denotes the 
 *           index of the unknown quantity to set (matrix row), J denotes
 *           the index of the distribution function to evaluate, and V is
 *           a scalar value to weight the distribution function with.
 */
void PXiExternalLoss::__SetElements(
    std::function<void(const len_t, const len_t, const real_t)> f
) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    const real_t *VpVol = this->grid->GetVpVol();
    
    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *mg = this->distributionGrid->GetMomentumGrid(ir);
        const len_t
            np  = mg->GetNp1(),
            nxi = mg->GetNp2();

        const real_t
            *dp    = mg->GetDp1(),
            *dp_f  = mg->GetDp1_f(),
            *dxi   = mg->GetDp2(),
            *Vp    = this->distributionGrid->GetVp(ir),
            *Vp_fp = this->distributionGrid->GetVp_f1(ir);

        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        real_t dd = 0;
        if (this->boundaryCondition == BC_DPHI_CONST)
            dd = dp[np-1] / dp[np-2];

        for (len_t j = 0; j < nxi; j++) {
            // Select correct indices/volume elements, depending on
            // whether we're building
            len_t idx1, idx2 = offset + j*np + (np-1);
            real_t iVd;
            if (this->boundary == BOUNDARY_FLUID) {
                idx1 = ir;
                // The fluxes should be reversed on this side of the boundary,
                // and we implement this reversal by negating this quantity
                // (as it is multiplied with the flux everywhere)
                iVd = -VpVol[ir] / dxi[j];
            } else {
                idx1 = idx2;
                iVd  = Vp[idx2-offset] * dp[np-1];
            }
            if(!iVd)
                continue;
            
            // Contribution from advection and PP diffusion
            if (this->boundaryCondition == BC_F_0) {
                real_t Vd = Vp_fp[j*(np+1) + np] / iVd;

                // the interpolation on the outermost cell interface is set to 
                // UPWIND: zero flux if negative advection, but free flow if positive. 
                real_t delta1 = (Ap[j*(np+1) + np]>0); 
                // Phi_{N_p+1/2}  -- f_{N_p+1} = 0
                f(idx1, idx2, delta1*Ap[j*(np+1) + np] * Vd);

                // Dpp
                f(idx1, idx2, Dpp[j*(np+1) + np]/dp_f[np-2] * Vd);

                if (j > 0 && j < nxi-1) {
                    f(idx1, idx2+np-1, -Dpx[j*(np+1) + np] / (2*dp_f[np-2]) * Vd);
                    f(idx1, idx2-np-1, +Dpx[j*(np+1) + np] / (2*dp_f[np-2]) * Vd);
                }
            } else {
                const real_t *delta1_0 = equation->GetInterpolationCoeff1(ir,np-1,j);
                const real_t *delta1_1 = equation->GetInterpolationCoeff1(ir,np-2,j);
                real_t Vd = Vp_fp[j*(np+1) + np-1] / iVd;

                // Phi_{N_p+1/2} = Phi_{N_p-1/2} + dd*(Phi_{N_p-1/2} - Phi_{N_p-3/2})
                // Phi_{N_p-1/2}
                for(len_t k=0; k<3; k++)
                    f(idx1, idx2+k-2, (1+dd)*delta1_0[k]*Ap[j*(np+1) + np-1] * Vd);
                
                // Phi_{N_p-3/2}
                for(len_t k=0; k<4; k++)
                    f(idx1, idx2+k-3, -dd*delta1_1[k]*Ap[j*(np+1) + np-1] * Vd);

                // Dpp
                f(idx1, idx2,   -(1+dd)*Dpp[j*(np+1) + np-1]/dp_f[np-2] * Vd);
                f(idx1, idx2-1,  (1+dd)*Dpp[j*(np+1) + np-1]/dp_f[np-2] * Vd);

                f(idx1, idx2-1,  dd*Dpp[j*(np+1) + np-2]/dp_f[np-3] * Vd);
                f(idx1, idx2-2, -dd*Dpp[j*(np+1) + np-2]/dp_f[np-3] * Vd);

                // Dpx
                if (j > 0 && j < nxi-1) {
                    f(idx1, idx2+np,   -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);
                    f(idx1, idx2+np-1, -(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);
                    f(idx1, idx2-np,   +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);
                    f(idx1, idx2-np-1, +(1+dd)*Dpx[j*(np+1) + np-1] / (dp_f[np-2]+dp_f[np-3]) * Vd);

                    f(idx1, idx2+np-1, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
                    f(idx1, idx2+np-2, +dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
                    f(idx1, idx2-np-1, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
                    f(idx1, idx2-np-2, -dd*Dpx[j*(np+1) + np-2] / (dp_f[np-3]+dp_f[np-4]) * Vd);
                }
            }
        }

        offset += np*nxi;
    }
}

