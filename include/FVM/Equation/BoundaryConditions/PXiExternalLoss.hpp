#ifndef _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_LOSS_HPP
#define _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_LOSS_HPP

#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class PXiExternalLoss : public BoundaryCondition {
    private:
        const Equation *equation;

    public:
        PXiExternalLoss(DREAM::FVM::Grid*, const DREAM::FVM::Equation*);
        virtual ~PXiExternalLoss();

        virtual bool Rebuild(const real_t, UnknownQuantityHandler*) override;

        virtual void AddToJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*) override;
        virtual void AddToMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;

        // Not implemented (not used)
        virtual void SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*) {}
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) {}
        virtual void SetVectorElements(real_t*, const real_t*) {}
    };
}

#endif/*_DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_LOSS_HPP*/
