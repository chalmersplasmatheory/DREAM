#ifndef _DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/Equations/AnalyticDistributionHottail.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include <gsl/gsl_roots.h>

namespace DREAM {
    class HottailRateTerm : public FVM::EquationTerm, public RunawaySourceTerm {
    private:
        void DeallocateAll();
        
    protected:
        AnalyticDistributionHottail *distHT;
        FVM::UnknownQuantityHandler *unknowns;
        real_t scaleFactor;
        real_t *pCrit = nullptr;
        real_t *gamma = nullptr;
        real_t *dGammaDPc = nullptr;
    public:
        HottailRateTerm(
            FVM::Grid*, AnalyticDistributionHottail*, 
            FVM::UnknownQuantityHandler*, real_t sf
        );
        ~HottailRateTerm();
        
        const real_t* GetRunawayRate() const { return gamma; }
        const real_t GetRunawayRate(const len_t ir) const { return gamma[ir]; }

        const real_t* GetHottailCriticalMomentum() const { return pCrit; }
        const real_t GetHottailCriticalMomentum(const len_t ir) const { return pCrit[ir]; }

        virtual bool GridRebuilt() override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}

        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override {return false;}
    };
}

#endif/*_DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_HPP*/
