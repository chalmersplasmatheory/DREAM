#ifndef _DREAM_EQUATIONS_TRANSPORT_PRESCRIBED_BC_HPP
#define _DREAM_EQUATIONS_TRANSPORT_PRESCRIBED_BC_HPP

#include <functional>
#include "DREAM/Equations/TransportPrescribed.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Equation/BoundaryCondition.hpp"

namespace DREAM {
    template<typename T>
    class TransportPrescribedBC : public FVM::BC::BoundaryCondition {
    private:
        TransportPrescribed<T> *transportOperator;

        void __SetElements(std::function<void(const len_t, const len_t, const real_t)>);
        real_t __GetSingleElement(const real_t, const real_t, const real_t);

    public:
        TransportPrescribedBC<T>(
            FVM::Grid*, TransportPrescribed<T>*
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
    real_t TransportPrescribedBC<FVM::AdvectionTerm>::__GetSingleElement(
        const real_t, const real_t, const real_t
    );
    template<>
    real_t TransportPrescribedBC<FVM::DiffusionTerm>::__GetSingleElement(
        const real_t, const real_t, const real_t
    );

    typedef TransportPrescribedBC<FVM::AdvectionTerm> TransportPrescribedAdvectiveBC;
    typedef TransportPrescribedBC<FVM::DiffusionTerm> TransportPrescribedDiffusiveBC;
}

#include "TransportPrescribedBC.tcc"

#endif/*_DREAM_EQUATIONS_TRANSPORT_PRESCRIBED_BC_HPP*/
