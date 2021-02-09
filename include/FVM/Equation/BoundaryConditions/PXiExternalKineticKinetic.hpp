#ifndef _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_KINETIC_KINETIC_HPP
#define _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_KINETIC_KINETIC_HPP

#include "FVM/Equation/BoundaryConditions/PXiAdvectionDiffusionBoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class PXiExternalKineticKinetic : public PXiAdvectionDiffusionBoundaryCondition {
    public:
        enum condition_type {
            TYPE_LOWER,     // B.C. at p=p0
            TYPE_UPPER,     // B.C. at p=pmax
            TYPE_DENSITY    // Flow of particles into fluid quantity
        };
    private:
        Grid *lowerGrid, *upperGrid;
        len_t id_f_low, id_f_upp;
        enum condition_type type;

        const real_t *fLow, *fUpp;

        void __SetElements(
            std::function<void(const len_t, const len_t, const real_t)>,
            std::function<void(const len_t, const len_t, const real_t)>
        );
        void __SetElements(
            std::function<void(const len_t, const len_t, const real_t)>,
            std::function<void(const len_t, const len_t, const real_t)>,
            const real_t *const*, const real_t *const*, const real_t *const*
        );

    public:
        PXiExternalKineticKinetic(
            DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            const DREAM::FVM::Operator*, const len_t, const len_t,
            enum condition_type
        );
        virtual ~PXiExternalKineticKinetic();

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

#endif/*_DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_KINETIC_KINETIC_HPP*/
