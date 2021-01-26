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
        struct altPcParams {len_t ir; real_t ncold; real_t Eterm; real_t tau; real_t lnL; IonHandler *ionHandler; AnalyticDistributionHottail *dist; FVM::RadialGrid *rGrid; real_t *fPointer; real_t *dfdpPointer;};

        enum OptionConstants::eqterm_hottail_mode type;
        real_t scaleFactor;
        AnalyticDistributionHottail *distHT;
        FVM::UnknownQuantityHandler *unknowns;
        CoulombLogarithm *lnL;
        const len_t 
            id_ncold,
            id_Efield,
            id_tau;

        real_t *pcAlt_prev = nullptr;
        real_t *pcAlt = nullptr;
        real_t *gamma = nullptr;
        real_t *dGammaDPc = nullptr;
        gsl_root_fdfsolver *fdfsolver;
        gsl_function_fdf gsl_altPcFunc;
        altPcParams gsl_altPcParams;

        real_t tPrev = -1.0;
        real_t dt;
        const real_t RELTOL_FOR_PC = 1e-10;
        const real_t ABSTOL_FOR_PC = 0.0;
        

        static real_t altPcFunc(real_t p, void *par);
        static real_t altPcFunc_df(real_t p, void *par);
        static void altPcFunc_fdf(real_t p, void *par, real_t *f, real_t *df);

        real_t evaluateAltCriticalMomentum(len_t ir, real_t *f=nullptr, real_t *dfdp=nullptr);
        real_t evaluatePartialAltCriticalMomentum(len_t ir, len_t derivId);
        void DeallocateAll();
    public:
        HottailRateTerm(
            FVM::Grid*, AnalyticDistributionHottail*, FVM::UnknownQuantityHandler*,
            IonHandler*, CoulombLogarithm*,
            enum OptionConstants::eqterm_hottail_mode, real_t scaleFactor=1.0
        );
        ~HottailRateTerm();
        
        const real_t* GetRunawayRate() const { return gamma; }
        const real_t GetRunawayRate(const len_t ir) const { return gamma[ir]; }

        const real_t* GetHottailCriticalMomentum() const { return type==OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC_ALT_PC ? pcAlt :  /* todo */ nullptr; }
        const real_t GetHottailCriticalMomentum(const len_t ir) const { return type==OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC_ALT_PC ? pcAlt[ir] :  /* todo */ NAN; }


        virtual bool GridRebuilt() override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }   /* XXX TODO XXX */
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_HPP*/
