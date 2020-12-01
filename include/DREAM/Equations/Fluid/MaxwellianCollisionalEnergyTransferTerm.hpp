#ifndef _DREAM_EQUATION_FLUID_MAXWELLIAN_COLLISIONAL_ENERGY_TRANSFER_TERM_HPP
#define _DREAM_EQUATION_FLUID_MAXWELLIAN_COLLISIONAL_ENERGY_TRANSFER_TERM_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "DREAM/IonHandler.hpp"

namespace DREAM{
    class MaxwellianCollisionalEnergyTransferTerm : public FVM::EquationTerm {
    private:
        real_t constPreFactor;
        const len_t index_i, index_j;
        bool isIon_i, isIon_j, isEI;
        FVM::UnknownQuantityHandler *unknowns;
        CoulombLogarithm *lnLambda;
        IonHandler *ionHandler;
        len_t Zi, Zj;
        real_t mi, mj;
        len_t 
            id_ions,
            id_Tcold,
            id_ncold,id_Wcold,
            id_Ni, id_Wi;

        struct CollisionQuantity::collqty_settings *lnLambda_settings;
        void GetParametersForSpecies(len_t ir, len_t index, bool isIon, real_t &n, real_t &W, real_t &nZ2);
    public:
        MaxwellianCollisionalEnergyTransferTerm(
            FVM::Grid *g, 
            const len_t index_i, bool isIon_i,
            const len_t index_j, bool isIon_j, 
            FVM::UnknownQuantityHandler *u, CoulombLogarithm *lnL, IonHandler *ionHandler, real_t scaleFactor=1.0
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
