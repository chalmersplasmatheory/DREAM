#ifndef _DREAM_EQUATION_MOLECULE_CHARGE_EXCHANGE_RATE_REACTION_HPP
#define _DREAM_EQUATION_MOLECULE_CHARGE_EXCHANGE_RATE_REACTION_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "DREAM/MoleculeHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class MoleculeChargeExchangeRateReaction : public IonEquationTerm<FVM::EquationTerm> {
    protected:
        enum SetMode {MATRIX, JACOBIAN};
        ADAS *adas;
        FVM::UnknownQuantityHandler *unknowns;
        RateHandler *ratehandler;
        len_t id_ions, id_n_cold, id_n_tot, id_T_cold;
        bool addFluidIonization; // the full ADAS ionization rate is added in this equation term
        bool addFluidJacobian;   // only the jacobian of the ionization is set with this term
        real_t **Rate = nullptr;
        real_t **Rate_N = nullptr;
        real_t **Rate_T = nullptr;

        len_t FindRelevantMolecularReactions();
        std::vector<MolecularReaction> chargeExchangeReactions;
       
    public:
        MoleculeChargeExchangeRateReaction(
            FVM::Grid*, IonHandler*, const len_t,ADAS*, 
            FVM::UnknownQuantityHandler*, RateHandler*, bool,bool,bool
        );
        virtual ~MoleculeChargeExchangeRateReaction();

        void AllocateRateCoefficients(const len_t allocationSize);
        void DeallocateRateCoefficients();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 2 * chargeExchangeReactions.size(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
            {
                return GetNumberOfNonZerosPerRow() + 2;
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

    
    #endif /** _DREAM_EQUATION_MOLECULE_CHARGE_EXCHANGE_RATE_REACTION_HPP */