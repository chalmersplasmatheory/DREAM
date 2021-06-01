#ifndef _DREAM_EQUATION_FLUID_CURRENT_DENSITY_FROM_ANALYTIC_RE_HPP
#define _DREAM_EQUATION_FLUID_CURRENT_DENSITY_FROM_ANALYTIC_RE_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/AnalyticDistributionRE.hpp"
#include "DREAM/Equations/REPitchDistributionAveragedBACoeff.hpp"
#include <gsl/gsl_integration.h>
namespace DREAM {
    class CurrentDensityFromAnalyticRE : public FVM::EquationTerm {
    public:
        struct integrandParameters {
            len_t ir; 
            AnalyticDistributionRE *distRE; 
            REPitchDistributionAveragedBACoeff *AveragedXiTerm;
            len_t derivId;
            len_t nMultiple;
        };
    private:
        AnalyticDistributionRE *distRE;
        FVM::UnknownQuantityHandler *unknowns;
        const len_t 
            id_Efield,
            id_ncold;
        real_t scaleFactor;

        REPitchDistributionAveragedBACoeff *AveragedXiTerm;
        
        real_t *currentDensity=nullptr;
        real_t *partialCurrentDensity=nullptr;
        gsl_integration_workspace *gsl_ad;
        gsl_function gsl_func, gsl_func_partial;
        integrandParameters gsl_params;

        static real_t currentDensityIntegrand(real_t p, void*);
        static real_t partialCurrentDensityIntegrand(real_t p, void*);

        void setPartialCurrentDensity(len_t derivId, len_t nMultiple);
        void Deallocate();
    public:
        CurrentDensityFromAnalyticRE(
            FVM::Grid*, FVM::UnknownQuantityHandler *u, 
            AnalyticDistributionRE*, real_t scaleFactor=1.0
        );
        virtual ~CurrentDensityFromAnalyticRE();

        virtual len_t GetNumberOfNonZerosPerRow() const override
            { return 1; }

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_CURRENT_DENSITY_FROM_ANALYTIC_RE_HPP*/
