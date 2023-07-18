#ifndef _DREAM_EQUATION_FLUID_AVALANCHE_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_AVALANCHE_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

namespace DREAM {
    class AvalancheCurrentDensityFromAnalyticalDistributionFunction : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        const real_t scaleFactor;


        len_t id_Efield;




        struct integrandRosenbluthPutvinskiParams {

        };


        static real_t integrandHesslow(real_t, void*);
        // static real_t integrandRosenbluthPutvinski(real_t, void*);
        real_t evaluateAvalancheRunawaysMeanVelocity(len_t);


    protected:

        gsl_integration_workspace *gsl_ad_w = nullptr;
        int QAG_KEY = GSL_INTEG_GAUSS31;

        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;

    public:
        AvalancheCurrentDensityFromAnalyticalDistributionFunction(
            FVM::Grid*, FVM::UnknownQuantityHandler*, RunawayFluid*, real_t scaleFactor=1.0);
        ~AvalancheCurrentDensityFromAnalyticalDistributionFunction();
    };
}


#endif /*_DREAM_EQUATION_FLUID_AVALANCHE_CURRENT_DENSITY_FROM_ANALYTICAL_DISTRIBUTION_FUNCTION_HPP*/
