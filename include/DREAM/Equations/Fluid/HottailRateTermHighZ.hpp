#ifndef _DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_HIGH_Z_HPP
#define _DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_HIGH_Z_HPP

#include "DREAM/Equations/Fluid/HottailRateTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include <gsl/gsl_roots.h>

namespace DREAM {
    class HottailRateTermHighZ : public HottailRateTerm {
    private:
        struct PcParams {len_t ir; real_t ncold; real_t Eterm; real_t tau; real_t lnL; IonHandler *ionHandler; AnalyticDistributionHottail *dist; FVM::RadialGrid *rGrid; real_t F; real_t dFdp;};

        CoulombLogarithm *lnL;
        const len_t 
            id_ncold,
            id_Efield,
            id_tau;

        real_t *pCrit_prev = nullptr;

        gsl_root_fdfsolver *fdfsolver;
        gsl_function_fdf gsl_func;
        PcParams gsl_params;

        real_t tPrev = -1.0;
        real_t dt;
        const real_t RELTOL_FOR_PC = 1e-10;
        const real_t ABSTOL_FOR_PC = 0.0;
        

        static real_t PcFunc(real_t p, void *par);
        static real_t PcFunc_df(real_t p, void *par);
        static void PcFunc_fdf(real_t p, void *par, real_t *f, real_t *df);

        real_t evaluateCriticalMomentum(len_t ir, real_t &f, real_t &dfdp);
        real_t evaluatePartialCriticalMomentum(len_t ir, len_t derivId);
        void Deallocate();
    public:
        HottailRateTermHighZ(
            FVM::Grid*, AnalyticDistributionHottail*, FVM::UnknownQuantityHandler*,
            IonHandler*, CoulombLogarithm*, real_t scaleFactor=1.0
        );
        ~HottailRateTermHighZ();
        
        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        
        virtual const real_t *GetElectricField() override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_HIGH_Z_HPP*/
