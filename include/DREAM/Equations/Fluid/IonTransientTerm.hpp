#ifndef _DREAM_EQUATIONS_FLUID_ION_TRANSIENT_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_ION_TRANSIENT_TERM_HPP

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"

namespace DREAM {
    class IonTransientTerm : public IonEquationTerm<FVM::EquationTerm> {
    private:
        real_t dt;
        len_t unknownId;

        // Ion densities in the previous time step
        real_t *xn;

    public:
        IonTransientTerm(FVM::Grid*, IonHandler*, const len_t iIon, const len_t unknownId);

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual bool SetCSJacobianBlock(
            const len_t, const len_t, FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;

        virtual void SetCSMatrixElements(
            FVM::Matrix*, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
        
        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
    };
}

#endif/*_DREAM_EQUATIONS_FLUID_ION_TRANSIENT_TERM_HPP*/
