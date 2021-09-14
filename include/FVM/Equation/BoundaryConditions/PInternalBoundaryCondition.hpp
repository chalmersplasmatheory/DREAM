#ifndef _DREAM_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP
#define _DREAM_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP

#include "FVM/config.h"
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class PInternalBoundaryCondition : public BoundaryCondition {
    protected:
        bool allZero = false;
        real_t **VpS = nullptr;
        len_t nr, *nxi;

        void _AddToVector(real_t*);

    public:
        PInternalBoundaryCondition(Grid*);

        void AllocateFluxes();
        void DeallocateFluxes();
        real_t& Flux(const len_t ir, const len_t j) { return this->VpS[ir][j]; }
        
        virtual bool GridRebuilt() override;
        virtual bool Rebuild(const real_t, UnknownQuantityHandler*) override;

        virtual bool AddToJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override {return false;}
        virtual void AddToMatrixElements(Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override {return false;}
        virtual void SetMatrixElements(Matrix*, real_t*) override {};
        virtual void SetVectorElements(real_t*, const real_t*) override {}
    };
}

#endif/*_DREAM_FVM_BC_P_INTERNAL_BOUNDARY_CONDITION_HPP*/
