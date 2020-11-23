#ifndef _DREAM_EQUATIONS_TRANSPORT_BC_HPP
#define _DREAM_EQUATIONS_TRANSPORT_BC_HPP

#include <functional>
#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    template<typename T>
    class TransportBC : public FVM::BC::BoundaryCondition {
    private:
        T *transportOperator;

        void __SetElements(std::function<void(const len_t, const len_t, const real_t)>);
        real_t __GetSingleElement(const real_t, const real_t, const real_t);
        const real_t *GetCoefficient(const len_t);

    public:
        TransportBC<T>(
            FVM::Grid*, T*
        );

        // Rebuilding is handled by the
        virtual bool Rebuild(const real_t, FVM::UnknownQuantityHandler*) override { return false; }

        virtual void AddToJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override {}
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}
    };

    template<>
    real_t TransportBC<FVM::AdvectionTerm>::__GetSingleElement(
        const real_t, const real_t, const real_t
    );
    template<>
    real_t TransportBC<FVM::DiffusionTerm>::__GetSingleElement(
        const real_t, const real_t, const real_t
    );

    template<>
    const real_t *TransportBC<FVM::AdvectionTerm>::GetCoefficient(const len_t);
    template<>
    const real_t *TransportBC<FVM::DiffusionTerm>::GetCoefficient(const len_t);

    typedef TransportBC<FVM::AdvectionTerm> TransportAdvectiveBC;
    typedef TransportBC<FVM::DiffusionTerm> TransportDiffusiveBC;
}

#include "TransportBC.tcc"

#endif/*_DREAM_EQUATIONS_TRANSPORT_BC_HPP*/
