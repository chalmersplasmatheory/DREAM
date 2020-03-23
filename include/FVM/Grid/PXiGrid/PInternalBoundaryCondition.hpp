#ifndef _DREAM_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP
#define _DREAM_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP

#include "FVM/config.h"
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM::FVM::BC {
    class PInternalBoundaryCondition : public BoundaryCondition {
    private:
        bool allZero = false;
        real_t **p2S = nullptr;
        len_t nr, *nxi;

    public:
        PInternalBoundaryCondition(RadialGrid *rg) : BoundaryCondition(rg) {};

        void AllocateFluxes();
        void DeallocateFluxes();
        real_t& Flux(const len_t ir, const len_t j) { return this->p2S[ir][j]; }
        
        virtual bool GridRebuilt() override;
        virtual bool Rebuild(const real_t) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
    };
}

#endif/*_DREAM_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP*/
