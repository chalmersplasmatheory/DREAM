#ifndef _DREAM_REACTION_ENERGY_TRANSFER_TERM_HPP
#define _DREAM_REACTION_ENERGY_TRANSFER_TERM_HPP

#include <vector>
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/Fluid/RateHandler.hpp"

namespace DREAM {
    class ReactionEnergyTransferTerm : public FVM::EquationTerm {
    private:
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;
        RateHandler *ratehandler;
        len_t iz, id_Te;
        bool neutralEquation;   // true if this term is in the W_n equation, false if in W_i
        bool separateNeutrals;

        std::vector<len_t> reactions;
        std::vector<real_t> share;
        real_t iIonMass;

        real_t *Energy = nullptr;

    public:
        ReactionEnergyTransferTerm(
            FVM::Grid*, len_t, FVM::UnknownQuantityHandler*, IonHandler*, RateHandler*,
            bool neutralEquation, bool includeFluidREIonization
        );
        virtual ~ReactionEnergyTransferTerm();

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override {
            return 2*reactions.size() + ions->GetNZ() + 2;
        }
    };
}

#endif/*_DREAM_REACTION_ENERGY_TRANSFER_TERM_HPP*/
