#ifndef _DREAM_FVM_BOUNDARY_CONDITION_HPP
#define _DREAM_FVM_BOUNDARY_CONDITION_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class BoundaryCondition {
    protected:
        Grid *grid;

    public:
        BoundaryCondition(Grid *g) : grid(g) {};
        virtual ~BoundaryCondition() {}

        virtual bool GridRebuilt() { return false; }

        virtual bool Rebuild(const real_t t, UnknownQuantityHandler*) = 0;

        virtual void AddToJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) = 0;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) = 0;
        virtual void AddToMatrixElements(Matrix*, real_t*) =0;
        virtual void SetMatrixElements(Matrix*, real_t*) = 0;
        virtual void AddToVectorElements(real_t*, const real_t*) = 0;
        virtual void SetVectorElements(real_t*, const real_t*) = 0;
    };
}

#endif/*_DREAM_FVM_BOUNDARY_CONDITION_HPP*/
