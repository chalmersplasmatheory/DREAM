#ifndef _DREAM_EQUATION_FLUID_MAXWELLIAN_COLLISIONAL_ENERGY_TRANSFER_TERM_HPP
#define _DREAM_EQUATION_FLUID_MAXWELLIAN_COLLISIONAL_ENERGY_TRANSFER_TERM_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"

namespace DREAM{
    class MaxwellianCollisionalEnergyTransferTerm : public FVM::EquationTerm {
    private:
        real_t constPreFactor;
        const len_t 
            id_ni, id_nj, 
            id_Wi, id_Wj,
            Zi, Zj,
            offset_i, offset_j;
        const real_t mi, mj;
        FVM::UnknownQuantityHandler *unknowns;
        CoulombLogarithm *lnLambda;
        bool isEI;
        len_t id_ions;
        len_t id_Tcold;
        struct CollisionQuantity::collqty_settings *lnLambda_settings;

    public:
        MaxwellianCollisionalEnergyTransferTerm(
            FVM::Grid*, 
            const len_t id_ni, const len_t id_Wi, const len_t Zi, const real_t mi, const len_t offset_i, 
            const len_t id_nj, const len_t id_Wj, const len_t Zj, const real_t mj, const len_t offset_j, 
            FVM::UnknownQuantityHandler*, CoulombLogarithm*, bool isEI, real_t scaleFactor=1.0
        );
        ~MaxwellianCollisionalEnergyTransferTerm() {delete lnLambda_settings;}
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 5 + isEI?1:unknowns->GetUnknown(id_ions)->NumberOfMultiples(); }
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) {}

        virtual void SetJacobianBlock(const len_t uqtyI, const len_t derivId, FVM::Matrix*, const real_t*);
        virtual void SetMatrixElements(FVM::Matrix*, real_t*);
        virtual void SetVectorElements(real_t*, const real_t*);
    };
}

#endif/*_DREAM_EQUATION_FLUID_MAXWELLIAN_COLLISIONAL_ENERGY_TRANSFER_TERM_HPP*/
