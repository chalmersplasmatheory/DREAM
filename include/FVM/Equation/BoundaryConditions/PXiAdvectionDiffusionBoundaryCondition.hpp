#ifndef _DREAM_FVM_BOUNDARY_CONDITION_P_XI_ADVECTION_DIFFUSION_HPP
#define _DREAM_FVM_BOUNDARY_CONDITION_P_XI_ADVECTION_DIFFUSION_HPP

#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"


namespace DREAM::FVM::BC {
    class PXiAdvectionDiffusionBoundaryCondition : public BoundaryCondition {
    public:
        // this is invoked when setting off-diagonal jacobian 
        // contributions on the radial flux grid 
        enum jacobian_interp_mode {
            NO_JACOBIAN         = 1, // used to SetVector and SetMatrixElements
            JACOBIAN_SET_LOWER  = 2, // sets offset contribution for radial flux  
            JACOBIAN_SET_CENTER = 3, // sets diagonal jacobian contributions
            JACOBIAN_SET_UPPER  = 4  // sets offset contribution for radial flux
        };

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
            const real_t *const*, const real_t *const* df1,
            const real_t *const* df2, const real_t *const* ddrr,
            const real_t *const* dd11, const real_t *const* dd12,
            const real_t *const* dd21, const real_t *const* dd22,
            jacobian_interp_mode set_mode=NO_JACOBIAN
        ) = 0;

        void AddPartialJacobianContributions(const len_t, const len_t, Matrix*, const real_t*, bool);
        void SetPartialJacobianContribution(
            const int_t, const len_t, Matrix*, const real_t*,
            jacobian_interp_mode,
            const real_t *const* dfr, const real_t *const* df1,
            const real_t *const* df2, const real_t *const* ddrr,
            const real_t *const* dd11, const real_t *const* dd12,
            const real_t *const* dd21, const real_t *const* dd22
        );
        void ResetJacobianColumn();
    };
}

#endif/*_DREAM_FVM_BOUNDARY_CONDITION_P_XI_ADVECTION_DIFFUSION_HPP*/
