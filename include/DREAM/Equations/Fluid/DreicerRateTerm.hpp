#ifndef _DREAM_EQUATION_FLUID_DREICER_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_DREICER_RATE_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class DreicerRateTerm : public FVM::EquationTerm {
    public:
        enum dreicer_type {
            CONNOR_HASTIE,  // Connor-Hastie runaway rate [Connor & Hastie, NF 15 (1975)]
            NEURAL_NETWORK  // Hesslow's neural network [Hesslow et al, JPP 85 (2019)]
        };

    private:
        RunawayFluid *REFluid;
        IonHandler *ions;
        len_t id_n_re;
        enum dreicer_type type = CONNOR_HASTIE;
        real_t scaleFactor=1.0;

        len_t id_E_field, id_n_cold;

        // Runaway rate
        real_t *gamma;

    public:
        DreicerRateTerm(
            FVM::Grid*, FVM::UnknownQuantityHandler*,
            RunawayFluid*, IonHandler*, enum dreicer_type,
            real_t scaleFactor=1.0
        );
        ~DreicerRateTerm();

        void AllocateGamma();
        void DeallocateGamma();
        
        real_t ConnorHastie(const len_t, const real_t, const real_t, const real_t);

        virtual bool GridRebuilt() override;
        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 1; }   /* XXX TODO XXX */
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_DREICER_RATE_TERM_HPP*/
