#ifndef _DREAM_ION_BOUNDARY_CONDITION_HPP
#define _DREAM_ION_BOUNDARY_CONDITION_HPP

#include <utility>
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM {
    template<class T>
    class IonBoundaryCondition : public T {
    protected:
        IonHandler *ions;

        // Index of ion species considered and its charge
        len_t iIon, Zion;

    public:
        IonBoundaryCondition(
            FVM::Grid*, IonHandler*, const len_t iIon,
            Args&& ... args
        );
        virtual ~IonBoundaryCondition();

        /**
         * Implemented in BoundaryCondition and empty here:
         *
         * virtual len_t GetNumberOfNonZerosPerRow() const = 0;
         * virtual len_t GetNumberOfNonZerosPerRow_jac() const = 0;
         * virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) = 0;
         */

        // Replace these 'Add/Set...' methods...
        virtual bool AddToJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        // ...with these 'AddCS/SetCS...' (CS = Charge State) methods.
        virtual bool AddCSToJacobianBlock(
            const len_t, const len_t, FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) = 0;
        virtual bool SetCSJacobianBlock(
            const len_t, const len_t, FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) = 0;
        virtual void AddCSToMatrixElements(
            FVM::Matrix*, real_t*, const len_t iIon,
            const len_t Z0, const len_t rOffset
        ) = 0;
        virtual void SetCSToMatrixElements(
            FVM::Matrix*, real_t*, const len_t iIon,
            const len_t Z0, const len_t rOffset
        ) = 0;
        virtual void AddCSToVectorElements(
            real_t*, const real_t*, const len_t iIon,
            const len_t Z0, const len_t rOffset
        ) = 0;
        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon,
            const len_t Z0, const len_t rOffset
        ) = 0;
    };

    #include "IonBoundaryCondition.tcc"
}

#endif/*_DREAM_ION_BOUNDARY_CONDITION_HPP*/
