#ifndef _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_LOSS_HPP
#define _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_LOSS_HPP

#include <functional>
#include "FVM/Equation/BoundaryConditions/PXiAdvectionDiffusionBoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class PXiExternalLoss : public PXiAdvectionDiffusionBoundaryCondition {
    public:
        enum bc_type {
            // Set f=0 at p=p_{Np+1} (the point after pmax,
            // just outside the computational domain)
            //
            // (This B.C. can help stabilize the system in
            // situations where the runaway region is > pmax)
            BC_F_0=1,

            // Set Phi_{N_p+1/2} = Phi_{N_p-1/2} i.e. that the
            // momentum-space flux is constant between the two
            // last p points of the momentum-space grid.
            // 
            // This B.C. is generally stable as long as the
            // runaway region is within the computational domain
            // (i.e. p_c < pmax)
            BC_PHI_CONST=2,

            // Set d Phi_{N_p+1/2} / dp = d Phi_{N_p-1/2} / dp.
            //
            // This B.C. makes physically some more sense than the
            // others, but is known to cause instabilities near
            // p=pmax. (Note that the there is no optimal way of
            // doing the B.C., and so this comment should NOT
            // be read to imply that the other B.C.s are unphysical;
            // they are just as physical as this condition.
            BC_DPHI_CONST=3
        };

        // Specifies which side of the boundary the condition is
        // applied to. Either to the fluid side (where we should
        // accumulate the influx in a fluid quantity), or the kinetic
        // side (where we must calculate the outflux of each cell
        // separately)
        enum boundary_type {
            // B.C. applied to fluid quantity
            BOUNDARY_FLUID,
            // B.C. applied to kinetic quantity
            BOUNDARY_KINETIC
        };
    private:
        FVM::Grid *distributionGrid=nullptr;

        len_t fId;

        enum bc_type boundaryCondition = BC_PHI_CONST;
        enum boundary_type boundary    = BOUNDARY_KINETIC;

        void __SetElements(std::function<void(const len_t, const len_t, const real_t)>);
        void __SetElements(
            std::function<void(const len_t, const len_t, const real_t)>,
            const real_t *const*, const real_t *const*, const real_t *const*
        );

    public:
        PXiExternalLoss(
            DREAM::FVM::Grid*, const DREAM::FVM::Operator*, const len_t,
            DREAM::FVM::Grid *distributionGrid=nullptr,
            enum boundary_type=BOUNDARY_KINETIC, enum bc_type bc=BC_PHI_CONST
        );

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual bool Rebuild(const real_t, UnknownQuantityHandler*) override;

        virtual void AddToJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
        virtual void AddToVectorElements_c(
            real_t*, const real_t*,
            const real_t *const*, const real_t *const* df1,
            const real_t *const* df2, const real_t *const* ddrr,
            const real_t *const* dd11, const real_t *const* dd12,
            const real_t *const* dd21, const real_t *const* dd22,
            jacobian_interp_mode set_mode=NO_JACOBIAN
        ) override;

        // Not implemented (not used)
        virtual void SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override {}
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}
    };
}

#endif/*_DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_LOSS_HPP*/
