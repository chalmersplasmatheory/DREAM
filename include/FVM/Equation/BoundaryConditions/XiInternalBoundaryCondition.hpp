#ifndef _DREAM_FVM_BC_XI_BOUNDARY_CONDITION_HPP
#define _DREAM_FVM_BC_XI_BOUNDARY_CONDITION_HPP

#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class XiInternalBoundaryCondition : public BoundaryCondition {
    public:
        XiInternalBoundaryCondition(Grid *g) : BoundaryCondition(g) {}

        bool Rebuild(const real_t, UnknownQuantityHandler*) override;

        virtual bool AddToJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override {return false;}
        virtual void AddToMatrixElements(Matrix*, real_t*) override {}
        virtual void AddToVectorElements(real_t*, const real_t*) override {}
        virtual bool SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override {return false;}
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override {}
    };
}

#endif/*_DREAM_FVM_BC_XI_BOUNDARY_CONDITION_HPP*/
