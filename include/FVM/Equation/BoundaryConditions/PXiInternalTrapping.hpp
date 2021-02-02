#ifndef _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_INTERNAL_TRAPPING_HPP
#define _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_INTERNAL_TRAPPING_HPP

#include <functional>
#include "FVM/Equation/BoundaryConditions/PXiAdvectionDiffusionBoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM::BC {
    class PXiInternalTrapping : public PXiAdvectionDiffusionBoundaryCondition {
    private:
        // List of indices corresponding to trapped particles
        // with negative xi0, at each radius
        PetscInt *nTrappedNegXi_indices=nullptr;   // size nr
        PetscInt **trappedNegXi_indices=nullptr;   // size nr x nTrappedNegXi_indices[ir]
        // The indices of the positive xi0 corresponding
        // (approximately) to the xi0 in 'trappedNegXi_indices'...
        PetscInt **trappedPosXi_indices=nullptr;   // size nr x nTrappedNegXi_indices[ir]

        // List of indices where radial fluxes should be mirrored
        PetscInt *nTrappedNegXiRadial_indices=nullptr;   // size nr
        PetscInt **trappedNegXiRadial_indices=nullptr;   // size nr x nTrappedNegXi_indices[ir]
        // The indices of the cells with positive xi0 containing
        // -xi0 from 'trappedNegXiRadial_indices'...
        PetscInt **trappedPosXiRadial_indices=nullptr;   // size nr x nTrappedNegXi_indices[ir]

        
        // Total number of rows in Matrix and Jacobian that should be reset
        len_t nRowsToReset;
        // Indices of all rows that should be reset in Matrix and Jacobian
        PetscInt *rowsToReset = nullptr;
        void _addElements(
            std::function<void(const len_t, const len_t, const real_t)>,
            const real_t *const*, const real_t *const*, const real_t *const*,
            const real_t *const*, const real_t *const*
        );
        len_t _setElements(
            const len_t, const len_t, std::function<void(const len_t, const len_t, const real_t)>
        );

        const real_t realeps = std::numeric_limits<real_t>::epsilon();
    public:
        PXiInternalTrapping(
            DREAM::FVM::Grid*, DREAM::FVM::Operator*
        );
        virtual ~PXiInternalTrapping();

        virtual bool GridRebuilt() override;
        virtual bool Rebuild(const real_t, UnknownQuantityHandler*) override;
        
        void DeallocateTrappedIndices();
        void LocateTrappedRegion();

        virtual void AddToJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(Matrix*, real_t*) override;
        virtual void AddToVectorElements_c(
            real_t*, const real_t*,
            const real_t *const* dfr, const real_t *const* df1=nullptr,
            const real_t *const* df2=nullptr, const real_t *const* ddrr=nullptr,
            const real_t *const* dd11=nullptr, const real_t *const* dd12=nullptr,
            const real_t *const* dd21=nullptr, const real_t *const* dd22=nullptr
        ) override;

        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_INTERNAL_TRAPPING_HPP*/
