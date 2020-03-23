#ifndef _DREAM_FVM_BOUNDARY_CONDITION_HPP
#define _DREAM_FVM_BOUNDARY_CONDITION_HPP

#include "FVM/config.h"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM::FVM::BC {
    class BoundaryCondition {
    protected:
        RadialGrid *grid;

    public:
        BoundaryCondition(RadialGrid *g) : grid(g) {};

        virtual bool GridRebuilt() { return false; }

        virtual bool Rebuild(const real_t t) = 0;
        virtual void SetMatrixElements(Matrix*) = 0;
    };
}

#endif/*_DREAM_FVM_BOUNDARY_CONDITION_HPP*/
