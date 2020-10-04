#ifndef _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_INTERNAL_TRAPPING_HPP
#define _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_INTERNAL_TRAPPING_HPP

#include <functional>
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM::BC {
    class PXiInternalTrapping : public BoundaryCondition {
    private:
        // List of indices corresponding to trapped particles
        // with negative xi0, at each radius
        PetscInt *nTrappedNegXi_indices=nullptr;   // size nr
        PetscInt **trappedNegXi_indices=nullptr;   // size nr x nTrappedNegXi_indices[ir]
        // The indices of the positive xi0 corresponding
        // (approximately) to the xi0 in 'trappedNegXi_indices'...
        PetscInt **trappedPosXi_indices=nullptr;   // size nr x nTrappedNegXi_indices[ir]

        len_t _setElements(
            const len_t, const len_t, std::function<void(const len_t, const len_t, const real_t)>
        );
    
    public:
        PXiInternalTrapping(
            DREAM::FVM::Grid*
        );
        virtual ~PXiInternalTrapping();

        virtual bool GridRebuilt() override;
        virtual bool Rebuild(const real_t, UnknownQuantityHandler*) override;
        
        void DeallocateTrappedIndices();
        void LocateTrappedRegion();

        virtual void AddToJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_INTERNAL_TRAPPING_HPP*/
