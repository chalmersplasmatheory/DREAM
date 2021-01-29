#ifndef _DREAM_FVM_BOUNDARY_CONDITION_P_XI_ADVECTION_DIFFUSION_HPP
#define _DREAM_FVM_BOUNDARY_CONDITION_P_XI_ADVECTION_DIFFUSION_HPP

#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"


namespace DREAM::FVM::BC {
    class PXiAdvectionDiffusionBoundaryCondition : public BoundaryCondition {
    protected:
        real_t *jacobianColumn;
        const Operator *oprtr;

    public:
        PXiAdvectionDiffusionBoundaryCondition(Grid*, const Operator*);
        virtual ~PXiAdvectionDiffusionBoundaryCondition();

        virtual bool GridRebuilt() override;

        virtual void AddToVectorElements(real_t*, const real_t*) override;
        virtual void AddToVectorElements_c(
            real_t*, const real_t*,
            const real_t *const*, const real_t *const* df2=nullptr,
            const real_t *const* dd11=nullptr, const real_t *const* dd12=nullptr,
            const real_t *const* dd21=nullptr, const real_t *const* dd22=nullptr
        ) = 0;

        void AddPartialJacobianContributions(const len_t, const len_t, Matrix*, const real_t*);
        void SetPartialJacobianContribution(
            const len_t, Matrix*, const real_t*,
            const real_t *const* df1=nullptr, const real_t *const* df2=nullptr,
            const real_t *const* dd11=nullptr, const real_t *const* dd12=nullptr,
            const real_t *const* dd21=nullptr, const real_t *const* dd22=nullptr
        );
        void ResetJacobianColumn();
    };
}

#endif/*_DREAM_FVM_BOUNDARY_CONDITION_P_XI_ADVECTION_DIFFUSION_HPP*/
