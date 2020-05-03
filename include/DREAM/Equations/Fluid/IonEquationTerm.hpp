#ifndef _DREAM_ION_EQUATION_TERM_HPP
#define _DREAM_ION_EQUATION_TERM_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM {
    class IonEquationTerm : public FVM::EquationTerm {
    private:
        IonHandler *ions;

    public:
        IonEquationTerm(FVM::Grid*, IonHandler*);
        virtual ~IonEquationTerm();

        /**
         * Implemented in EquationTerm and empty here:
         *
         * virtual len_t GetNumberOfNonZerosPerRow() const = 0;
         * virtual len_t GetNumberOfNonZerosPerRow_jac() const = 0;
         * virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) = 0;
         */

        // Replace these 'Set...' methods with...
        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        // ..these 'SetCS...' (CS = Charge State) methods.
        virtual void SetCSJacobianBlock(
            const len_t, const len_t, FVM::Matrix*,
            const len_t iIon, const len_t rOffset
        ) = 0;
        virtual void SetCSMatrixElements(
            FVM::Matrix*, real_t*, const len_t iIon, const len_t rOffset
        ) = 0;
        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon, const len_t rOffset
        ) = 0;
    };
}

#endif/*_DREAM_ION_EQUATION_TERM_HPP*/
