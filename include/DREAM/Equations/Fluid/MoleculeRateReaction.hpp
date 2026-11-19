#ifndef _DREAM_EQUATION_MOLECULE_RATE_REACTION_HPP
#define _DREAM_EQUATION_MOLECULE_RATE_REACTION_HPP

#include <vector>
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class MoleculeRateReaction : public IonEquationTerm<FVM::EquationTerm> {
    protected:
        FVM::UnknownQuantityHandler *unknowns;
        RateHandler *ratehandler;
        len_t id_ions;
        std::string iIonName;
        std::vector<len_t> reactions;

    public:
        MoleculeRateReaction(
            FVM::Grid*, IonHandler*, const len_t,
            FVM::UnknownQuantityHandler*, RateHandler*,
            bool includeFluidREIonization);

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 2 * reactions.size(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return ions->GetNzs() + 2*ions->GetNZ() + 3;}

        // The rates are rebuilt by the RateHandler
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}

        virtual bool SetCSJacobianBlock(
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

#endif/*_DREAM_EQUATION_MOLECULE_RATE_REACTION_HPP*/
