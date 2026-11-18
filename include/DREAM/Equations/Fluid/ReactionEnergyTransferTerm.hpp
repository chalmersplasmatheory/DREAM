#ifndef _DREAM_REACTION_ENERGY_TRANSFER_TERM_HPP
#define _DREAM_REACTION_ENERGY_TRANSFER_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/Fluid/RateHandler.hpp"

namespace DREAM {
    class ReactionEnergyTransferTerm : public FVM::EquationTerm {
    private:
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;
        RateHandler *ratehandler;
        len_t iz, id_ions, id_Wn, id_Nn, id_W, id_ne, id_Te, id_n_re;
        real_t mi;
        bool neutralEquation;   // true if this term is in the W_n equation, false if in W_i
        bool separateNeutrals;
        bool includeFluidREIonization;

        real_t getTotalMassReactionProducts(const MolecularReaction&) const;
        real_t getEnergyPerReaction(const MolecularReaction&, const len_t) const;

    public:
        ReactionEnergyTransferTerm(
            FVM::Grid*, len_t, FVM::UnknownQuantityHandler*, IonHandler*, RateHandler*,
            bool neutralEquation, bool includeFluidREIonization
        );
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override {
            // Per reaction: two reactant columns, two interpolation inputs, neutral
            // density and own heat. Plus the own heat for atomic ionization/recombination.
            return 6 * ratehandler->GetMolecularReactions().size() + 1;
        }
    };
}

#endif/*_DREAM_REACTION_ENERGY_TRANSFER_TERM_HPP*/
