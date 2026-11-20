#ifndef _DREAM_NEUTRAL_IONIZATION_ENERGY_TRANSFER_TERM_HPP
#define _DREAM_NEUTRAL_IONIZATION_ENERGY_TRANSFER_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/Fluid/RateHandler.hpp"

namespace DREAM {
    class NeutralIonizationEnergyTransferTerm : public FVM::EquationTerm {
    private:
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;
        RateHandler *ratehandler;
        len_t iz, id_ions, id_ne, id_Te, id_n_re, id_Wn, id_Nn, id_Wi, id_Ni;
        bool neutralEquation;   // true if this term is in the W_n equation, false if in W_i

        real_t reIonizationRate;

        real_t *Ion = nullptr, *Ion_N = nullptr, *Ion_T = nullptr;
        real_t *Rec = nullptr, *Rec_N = nullptr, *Rec_T = nullptr;

    public:
        NeutralIonizationEnergyTransferTerm(
            FVM::Grid*, len_t, FVM::UnknownQuantityHandler*, IonHandler*, RateHandler*,
            bool neutralEquation, bool includeFluidREIonization
        );
        virtual ~NeutralIonizationEnergyTransferTerm();

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 0; }
        // Jacobian columns per row: n_i (Z0=0 and Z0=1), n_cold, T_cold, n_re, W_n, N_n, W_i, N_i
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 9; }
    };
}

#endif/*_DREAM_NEUTRAL_IONIZATION_ENERGY_TRANSFER_TERM_HPP*/
