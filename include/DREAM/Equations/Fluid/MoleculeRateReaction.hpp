#ifndef _DREAM_EQUATION_MOLECULE_RATE_REACTION_HPP
#define _DREAM_EQUATION_MOLECULE_RATE_REACTION_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "DREAM/MoleculeHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class MoleculeRateReaction : public IonEquationTerm<FVM::EquationTerm> {
    protected:
        enum SetMode {MATRIX, JACOBIAN};
        ADAS *adas;
        FVM::UnknownQuantityHandler *unknowns;
        RateHandler *ratehandler;
        bool separateNeutrals;
        len_t id_ions, id_n_cold, id_n_re, id_n_tot, id_T_cold, id_Wi, id_Ni, id_Wn;
        bool includeFluidREIonization; // include IONIZATION_RUNAWAY reactions (ionization = FLUID_RE)
        real_t reionizScale;           // scale factor for runaway ionization (ions/reioniz_scale)
        real_t **Rate = nullptr;
        real_t **Rate_N = nullptr;
        real_t **Rate_T = nullptr;

        len_t FindRelevantMolecularReactions();
        std::vector<MolecularReaction> moleculeReactions;
        void GetReactantDensity(const MolecularReactionSpecies&, len_t &id, len_t &offset) const;

    public:
        MoleculeRateReaction(
            FVM::Grid*, IonHandler*, const len_t,
            FVM::UnknownQuantityHandler*, RateHandler*,
            bool includeFluidREIonization, real_t reionizScale);
        virtual ~MoleculeRateReaction();

        void AllocateRateCoefficients(const len_t allocationSize);
        void DeallocateRateCoefficients();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 2 * moleculeReactions.size(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
            {
                return ions->GetNzs() + 2*ions->GetNZ() + 3;
            }

        //virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

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


        //virtual void SetMatrixElement(const len_t, const len_t, const real_t) override;
        //virtual void SetJacobianElement(const len_t, const len_t, const real_t) override;
    };
}

    
    #endif /** _DREAM_EQUATION_MOLECULE_RATE_REACTION_HPP */