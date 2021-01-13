#ifndef _DREAM_EQUATION_FLUID_TRITIUM_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_TRITIUM_RATE_TERM_HPP

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class TritiumRateTerm : public IonEquationTerm<FVM::EquationTerm>, public RunawaySourceTerm {
    private:
        RunawayFluid *REFluid;
        real_t scaleFactor;

        void SetCSElements_internal(
            std::function<void(const len_t, const len_t, const real_t)>,
            const len_t, const len_t, const len_t
        );
    public:
        TritiumRateTerm(
            FVM::Grid*, IonHandler*, FVM::UnknownQuantityHandler*, const len_t,
            RunawayFluid*, real_t scaleFactor=1.0
        );

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
        
        virtual void SetCSJacobianBlock(
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

#endif /*_DREAM_EQUATION_FLUID_TRITIUM_RATE_TERM_HPP*/
