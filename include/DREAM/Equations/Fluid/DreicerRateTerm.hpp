#ifndef _DREAM_EQUATION_FLUID_DREICER_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_DREICER_RATE_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class DreicerRateTerm : public FVM::EquationTerm, public RunawaySourceTerm {
    public:
        enum dreicer_type {
            CONNOR_HASTIE_NOCORR,   // Connor-Hastie runaway rate (without corrections)
            CONNOR_HASTIE,  // Connor-Hastie runaway rate [Connor & Hastie, NF 15 (1975)]
            NEURAL_NETWORK  // Hesslow's neural network [Hesslow et al, JPP 85 (2019)]
        };

    private:
        FVM::UnknownQuantityHandler *unknowns;
        RunawayFluid *REFluid;
        IonHandler *ions;
        enum dreicer_type type = CONNOR_HASTIE;
        real_t scaleFactor=1.0;

        len_t id_E_field, id_n_cold, id_n_tot, id_T_cold;

        // Runaway rate
        real_t *gamma;
        // Derivative of runaway rate w.r.t. E/E_D,
        // times E/E_D
        real_t *EED_dgamma_dEED;

        const real_t *data_E_field, *data_n_cold, *data_n_tot, *data_T_cold;

    public:
        DreicerRateTerm(
            FVM::Grid*, FVM::UnknownQuantityHandler*,
            RunawayFluid*, IonHandler*, enum dreicer_type,
            real_t scaleFactor=1.0
        );
        ~DreicerRateTerm();

        void AllocateGamma();
        void DeallocateGamma();
        
        virtual bool GridRebuilt() override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }   /* XXX TODO XXX */
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_DREICER_RATE_TERM_HPP*/
