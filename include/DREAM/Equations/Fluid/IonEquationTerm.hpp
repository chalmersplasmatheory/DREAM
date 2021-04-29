#ifndef _DREAM_ION_EQUATION_TERM_HPP
#define _DREAM_ION_EQUATION_TERM_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/MomentQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM {
    template<class T>
    class IonEquationTerm : public T {
    protected:
        IonHandler *ions;
        // Index of ion species to which this equation term should
        // be applied
        len_t iIon;
        len_t Zion;   // Charge number of ion 'iIon'

    public:
        IonEquationTerm(FVM::Grid*, IonHandler*, const len_t iIon);
        IonEquationTerm(
            FVM::Grid*, FVM::Grid*, const len_t momentId, const len_t fId, 
            FVM::UnknownQuantityHandler*, real_t pThreshold, FVM::MomentQuantity::pThresholdMode,
             IonHandler*, const len_t iIon
        );
        virtual ~IonEquationTerm();

        /**
         * Implemented in EquationTerm and empty here:
         *
         * virtual len_t GetNumberOfNonZerosPerRow() const = 0;
         * virtual len_t GetNumberOfNonZerosPerRow_jac() const = 0;
         * virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) = 0;
         */

        // Replace these 'Set...' methods with...
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        // ..these 'SetCS...' (CS = Charge State) methods.
        virtual bool SetCSJacobianBlock(
            const len_t, const len_t, FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) = 0;
        virtual void SetCSMatrixElements(
            FVM::Matrix*, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) = 0;
        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) = 0;
    };

    #include "IonEquationTerm.tcc"
}

#endif/*_DREAM_ION_EQUATION_TERM_HPP*/
