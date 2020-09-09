#ifndef _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_UPPER_HPP
#define _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_UPPER_HPP

#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class PXiExternalKineticUpper : public BoundaryCondition {
    private:
        Grid *lowerGrid, *upperGrid;
        const Operator *oprtr;

        len_t id_f_low, id_f_upp;
        const real_t *fLow, *fUpp;

        void __SetElements(
            std::function<void(const len_t, const len_t, const real_t)>,
            std::function<void(const len_t, const len_t, const real_t)>
        );

    public:
        PXiExternalKineticUpper(
            DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            const DREAM::FVM::Operator*, const len_t, const len_t
        );
        virtual ~PXiExternalKineticUpper();

        virtual bool Rebuild(const real_t, UnknownQuantityHandler*) override;

        virtual void AddToJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;

        // Not implemented (not used)
        virtual void SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override {}
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}
    };
}

#endif/*_DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_UPPER_HPP*/
