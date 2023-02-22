#ifndef _DREAM_BOOTSTRAP_EQUATION_TERM_HPP
#define _DREAM_BOOTSTRAP_EQUATION_TERM_HPP

#include "FVM/Matrix.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/IonHandler.hpp"

namespace DREAM {
    class BootstrapEquationTerm : public FVM::EquationTerm {

    private:
        real_t scaleFactor;

        real_t *deltaX = nullptr; // x_{k+1} - x_{k-1}

        len_t id_X; // set in subclass, one of {id_ncold, id_Tcold, id_Ni, id_Wi}.
        len_t nZ;
        len_t nzs;

        len_t
            id_ncold,
            id_ions,
            id_Ni,
            id_Tcold,
            id_Wi;

    public:
        BootstrapEquationTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*, len_t);
        ~BootstrapEquationTerm();

        void AllocateDeltaX();
        void DeallocateDeltaX();

	virtual len_t GetNumberOfNonZerosPerRow() const { return 3; } // is this correct?

        // helpers
        void SetJacobianElement(len_t, FVM::Matrix*, len_t, len_t, real_t);
        void SetMatrixElement(FVM::Matrix, len_t, real_t);

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) 
override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        // these are required in any subclass
        virtual real_t GetCoefficient(len_t, len_t) = 0;
        virtual real_t GetPartialCoefficient(len_t, len_t, len_t, len_t) = 0;
    };
}


#endif /*_DREAM_BOOTSTRAP_EQUATION_TERM_HPP*/
