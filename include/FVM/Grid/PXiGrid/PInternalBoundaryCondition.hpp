#ifndef _TQS_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP
#define _TQS_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP

#include "FVM/Equation/BoundaryCondition.hpp"

namespace TQS::FVM::BC {
    class PInternalBoundaryCondition : public BoundaryCondition {
    private:
        bool allZero = false;
        real_t **p2S = nullptr;
        len_t nr, *nxi;

    public:
        PInternalBoundaryCondition(RadialGrid *rg) : BoundaryCondition(rg);

        real_t& Flux(const len_t ir, const len_t j) { return this->p2S[ir][j]; }
        
        virtual bool GridRebuilt() override;
        virtual bool Rebuild(const real_t) override;
        virtual void SetMatrixElements(Matrix*) override;
    };
}

#endif/*_TQS_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP*/
