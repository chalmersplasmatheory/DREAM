#ifndef _TQS_FVM_BC_XI_BOUNDARY_CONDITION_HPP
#define _TQS_FVM_BC_XI_BOUNDARY_CONDITION_HPP

#include "FVM/Equation/BoundaryCondition.hpp"

namespace TQS::FVM::BC {
    class XiInternalBoundaryCondition : public BoundaryCondition {
    public:
        XiInternalBoundaryCondition(RadialGrid *rg) : BoundaryCondition(rg);


    };
}

#endif/*_TQS_FVM_BC_XI_BOUNDARY_CONDITION_HPP*/
