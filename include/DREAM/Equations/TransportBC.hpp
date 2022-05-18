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
    public:
        enum bctype {
            TRANSPORT_BC_F0,            // Assume f=0 at r>rmax
            TRANSPORT_BC_DF_CONST       // Assume d^2 f / dr^2 = 0 at r>rmax (i.e. df/dr = const)
        };

        enum jacobian_interp_mode {
            NO_JACOBIAN         = 1,    // Used to 'SetVectorElements()' and 'SetMatrixElements()'
            JACOBIAN_SET_LOWER  = 2,    // Sets offset contribution for radial flux
            JACOBIAN_SET_CENTER = 3,    // Sets diagonal jacobian contributions
            JACOBIAN_SET_UPPER  = 4     // Sets offset contribution for radial flux
        };
    private:
        T *transportOperator;

        real_t *jacobianColumn;

        enum bctype type = TRANSPORT_BC_F0;

        void __SetElements(std::function<void(const len_t, const len_t, const real_t)>, jacobian_interp_mode);
        void __SetElements(std::function<void(const len_t, const len_t, const real_t)>, const real_t*, jacobian_interp_mode);
        real_t __GetSingleElement(const real_t, const real_t, const real_t);
        const real_t *const* GetCoefficient();
        const real_t *GetCoefficient(const len_t);
        const real_t *const* GetDiffCoefficient();
        const real_t *GetDiffCoefficient(const len_t);
        void SetPartialTerm(const len_t, const len_t);
        void ResetJacobianColumn();

    public:
        TransportBC<T>(
            FVM::Grid*, T*, enum bctype type=TRANSPORT_BC_F0
        );
        ~TransportBC<T>();

        virtual bool GridRebuilt() override;

        // Rebuilding is handled by the
        virtual bool Rebuild(const real_t, FVM::UnknownQuantityHandler*) override { return false; }

        virtual bool AddToJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;
        void AddToVectorElements_c(real_t*, const real_t*, const real_t *const*, jacobian_interp_mode);

        virtual bool SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override {return false;}
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}

        void SetPartialJacobianContribution(
            const int_t, const len_t, DREAM::FVM::Matrix*, const real_t*,
            jacobian_interp_mode, const real_t *const*
        );
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
    const real_t *const* TransportBC<FVM::AdvectionTerm>::GetCoefficient();
    template<>
    const real_t *const* TransportBC<FVM::DiffusionTerm>::GetCoefficient();

    template<>
    const real_t *TransportBC<FVM::AdvectionTerm>::GetCoefficient(const len_t);
    template<>
    const real_t *TransportBC<FVM::DiffusionTerm>::GetCoefficient(const len_t);

    template<>
    const real_t *const* TransportBC<FVM::AdvectionTerm>::GetDiffCoefficient();
    template<>
    const real_t *const* TransportBC<FVM::DiffusionTerm>::GetDiffCoefficient();

    template<>
    const real_t *TransportBC<FVM::AdvectionTerm>::GetDiffCoefficient(const len_t);
    template<>
    const real_t *TransportBC<FVM::DiffusionTerm>::GetDiffCoefficient(const len_t);

    template<>
    void TransportBC<FVM::AdvectionTerm>::SetPartialTerm(const len_t, const len_t);
    template<>
    void TransportBC<FVM::DiffusionTerm>::SetPartialTerm(const len_t, const len_t);

    typedef TransportBC<FVM::AdvectionTerm> TransportAdvectiveBC;
    typedef TransportBC<FVM::DiffusionTerm> TransportDiffusiveBC;
}

#include "TransportBC.tcc"

#endif/*_DREAM_EQUATIONS_TRANSPORT_BC_HPP*/
