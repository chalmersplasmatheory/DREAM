#ifndef _DREAM_EXOTERMIC_ENERGY_TRANSFER_TERM_HPP
#define _DREAM_EXOTERMIC_ENERGY_TRANSFER_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/Fluid/RateHandler.hpp"

namespace DREAM {
    class ExotermicEnergyTransferTerm : public FVM::EquationTerm {
    private:
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;
        RateHandler *ratehandler;
        len_t iz, id_ions, id_Wn, id_ne, id_Te;
        real_t mi;

        real_t getExotermicEnergy(const MolecularReaction&) const;
        real_t getTotalMassReactionProducts(const MolecularReaction&) const;

    public:
        ExotermicEnergyTransferTerm(
            FVM::Grid*, len_t, FVM::UnknownQuantityHandler*, IonHandler*, RateHandler*
        );
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override {
            // Two reactant columns, two interpolation inputs, and neutral density.
            return 5 * ratehandler->GetMolecularReactions().size();
        }
    };
}

#endif/*_DREAM_EXOTERMIC_ENERGY_TRANSFER_TERM_HPP*/
