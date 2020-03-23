#ifndef _DREAM_FVM_BC_XI_BOUNDARY_CONDITION_HPP
#define _DREAM_FVM_BC_XI_BOUNDARY_CONDITION_HPP

#include "FVM/Equation/BoundaryCondition.hpp"

namespace DREAM::FVM::BC {
    class XiInternalBoundaryCondition : public BoundaryCondition {
    public:
        XiInternalBoundaryCondition(RadialGrid *rg) : BoundaryCondition(rg);

        bool Rebuild(const real_t) override;
        void SetMatrixElements(Matrix*) override;
    };
}

#endif/*_DREAM_FVM_BC_XI_BOUNDARY_CONDITION_HPP*/
