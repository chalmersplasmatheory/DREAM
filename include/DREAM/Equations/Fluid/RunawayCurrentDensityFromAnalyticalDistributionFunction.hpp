#ifndef _DREAM_EQUATION_FLUID_RUNAWAY_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_RUNAWAY_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

namespace DREAM {
    class RunawayCurrentDensityFromAnalyticalDistributionFunction : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        const real_t scaleFactor;

        const real_t constPreFactor = Constants::c * Constants::ec;

        len_t id_Efield;
        len_t id_ntot;
        len_t id_ncold;

        real_t *Efield;
        const real_t *FSA_B;

        struct integrandParams {
            real_t beta;
            real_t pCrit;
        };

        static real_t integrand(real_t, void*);
        void evaluateMeanSpeed();

        gsl_integration_workspace *gsl_w;
        int GSL_WORKSPACE_SIZE = 1000;
        int QAG_KEY = GSL_INTEG_GAUSS31;

        real_t *u_re = nullptr;
        real_t *beta = nullptr;
        real_t *dPCrit = nullptr;

    protected:
        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;

    public:
        RunawayCurrentDensityFromAnalyticalDistributionFunction(
            FVM::Grid*, FVM::UnknownQuantityHandler*, RunawayFluid*,real_t scaleFactor=1.0);
        ~RunawayCurrentDensityFromAnalyticalDistributionFunction();
    };
}


#endif /*_DREAM_EQUATION_FLUID_RUNAWAY_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP*/
