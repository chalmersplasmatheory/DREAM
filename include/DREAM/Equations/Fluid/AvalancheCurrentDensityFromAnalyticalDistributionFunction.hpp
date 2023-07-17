#ifndef _DREAM_EQUATION_FLUID_AVALANCHE_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_AVALANCHE_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

namespace DREAM {
    class AvalancheCurrentDensityFromAnalyticalDistributionFunction : public FVM::DiagonalComplexTerm {
    private:
        const real_t scaleFactor;

        len_t id_ncold;
        len_t id_ntot;
        len_t id_Efield;

        /**
        * Parameter struct containing integrand parameters which is passed to a GSL function.
        */
        struct integrandHesslowParams {
            real_t ncold;
            real_t ntot;
            real_t lnLambda;
            real_t pceff;
            real_t nuDnuS;
        };
        struct integrandRosenbluthPutvinskiParams {
            
        };


        void InitializeGSLWorkspace();
        void DeallocateGSL();

    protected:

        gsl_integration_workspace *gsl_ad_w = nullptr;
        int QAG_KEY = GSL_INTEG_GAUSS31;

        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;

    public:
        AvalancheCurrentDensityFromAnalyticalDistributionFunction(
            FVM::Grid*, FVM::UnknownQuantityHandler*, real_t scaleFactor=1.0);
        ~AvalancheCurrentDensityFromAnalyticalDistributionFunction();
    };
}


#endif /*_DREAM_EQUATION_FLUID_AVALANCHE_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP*/
