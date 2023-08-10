#ifndef _DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_LOW_Z_HPP
#define _DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_LOW_Z_HPP

#include "DREAM/Equations/Fluid/HottailRateTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include <gsl/gsl_integration.h>

namespace DREAM {
    class HottailRateTermLowZ : public HottailRateTerm {
    private:
        struct ParamStruct {len_t ir; const real_t* ncold; const real_t* sp_cond; const real_t* tau; real_t* lnL; const real_t* j0; const real_t* Zeff; real_t* Ec; AnalyticDistributionHottail *dist; FVM::RadialGrid *rGrid;}; 
        
        struct dGammadu_Params {len_t ir; len_t id_unknown; real_t tau; real_t lnL; real_t ncold; AnalyticDistributionHottail *dist; real_t lowlim; len_t n;};

        CoulombLogarithm *lnL;
        IonHandler *ionHandler;
        RunawayFluid *runawayFluid;
        FVM::RadialGrid *rGrid;
        const len_t 
            id_ncold,
            id_Efield,
            id_tau,
            id_Tcold,
            id_johm,
            id_ni;

        real_t *Epar_prev = nullptr;
        real_t *Epar = nullptr;
        real_t *integraloff0inGamma;
        
        ParamStruct *params;

        real_t tPrev = -1.0;
        real_t dt;
        const real_t RELTOL_FOR_INT = 1e-7;
        const real_t ABSTOL_FOR_INT = 1e-250;
        const len_t nGslIntervals = 1000;
        
        void Deallocate();
        
        static real_t evaluate_f0(real_t p, void *par);
        static real_t partialIntegrandForEpar(real_t p);
        static real_t totalIntegrandForEpar(real_t p, void *par);
        real_t integralEpar(struct ParamStruct * intparams);
        real_t evaluate_Epar(struct ParamStruct * intparams);
        real_t integralf0(real_t Epar_ir, struct ParamStruct * intparams);
        
        real_t dEcdu(struct dGammadu_Params * pars);
        real_t dSpConddu(struct dGammadu_Params * pars);
        static real_t df0dtau(real_t p, void *par);
        static real_t totalIntegrandFor_dBoxIntdu(real_t p, void *par);
        real_t dBoxIntdu(struct dGammadu_Params * pars);
        real_t Intdf0du(struct dGammadu_Params * pars);
        real_t dGammadu(len_t ir, len_t id_unknown, len_t n);
        
        gsl_integration_workspace * w;
        
    public:
        HottailRateTermLowZ(
            FVM::Grid*, AnalyticDistributionHottail*, FVM::UnknownQuantityHandler*,
            IonHandler*, CoulombLogarithm*, RunawayFluid*, real_t scaleFactor=1.0
        );
        ~HottailRateTermLowZ();
        
        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        
        virtual const real_t *GetElectricField() override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_HOTTAIL_RATE_TERM_LOW_Z_HPP*/
