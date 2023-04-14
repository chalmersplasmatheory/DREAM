#ifndef _DREAM_BOOTSTRAP_EQUATION_TERM_HPP
#define _DREAM_BOOTSTRAP_EQUATION_TERM_HPP

#include "FVM/Matrix.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/BootstrapCurrent.hpp"

namespace DREAM {
    class BootstrapEquationTerm : public FVM::EquationTerm {

    private:
        OptionConstants::eqterm_bootstrap_bc bc;
        real_t scaleFactor;

        real_t *deltaX = nullptr; // x_{k+1} - x_{k-1}

        len_t id_X; // set in subclass, one of {id_ncold, id_Tcold, id_Ni, id_Wi}.
        len_t nzs;

        bool isInitialized;



        // helpers
        void rebuildBackwardsBC(FVM::UnknownQuantityHandler*);
        void rebuildZeroBC(FVM::UnknownQuantityHandler*);
        void setJacobianElement(len_t, FVM::Matrix*, len_t, len_t, real_t);
        void setJacobianElementBackwardsBC(len_t, FVM::Matrix*, len_t, len_t, real_t);
        void setJacobianElementZeroBC(len_t, FVM::Matrix*, len_t, len_t, real_t);
        void setMatrixElement(FVM::Matrix*, len_t, real_t);
        void setMatrixElementBackwardsBC(FVM::Matrix*, len_t, real_t);
        void setMatrixElementZeroBC(FVM::Matrix*, len_t, real_t);



    protected:
        BootstrapCurrent *bs;

        len_t nZ;

        len_t
            id_ncold,
            id_ions,
            id_Ni,
            id_Tcold,
            id_Wi;

    public:
        BootstrapEquationTerm(
            FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*,
            BootstrapCurrent*, OptionConstants::eqterm_bootstrap_bc, real_t
        );
        ~BootstrapEquationTerm();

        void AllocateDeltaX();
        void DeallocateDeltaX();

        void SetUnknownID(len_t id_X) { this->id_X = id_X; };

        // overrides of methods from EquationTerm
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual len_t GetNumberOfNonZerosPerRow() const { return 3; }; // is this correct?

        // these are required in any subclass of this class
        virtual real_t GetCoefficient(len_t, len_t) = 0;
        virtual real_t GetPartialCoefficient(len_t, len_t, len_t, len_t) = 0;
    };
}


#endif /*_DREAM_BOOTSTRAP_EQUATION_TERM_HPP*/
