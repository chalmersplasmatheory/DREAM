#ifndef _DREAMTESTS_DREAM_AVALANCHE_SOURCE_CH_HPP
#define _DREAMTESTS_DREAM_AVALANCHE_SOURCE_CH_HPP

#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
    class AvalancheSourceCH : public UnitTest {
    private:
        int QAG_KEY = GSL_INTEG_GAUSS61;
        gsl_integration_workspace *gsl_ws1, *gsl_ws2, *gsl_ws3;
        struct intFdistParams {real_t pin; real_t mag;};
        static real_t distributionFunction(real_t xi, void *par);
        struct intXiParams {real_t pin; real_t mag; int qag_key; gsl_integration_workspace *gslws3;};
        static real_t analyticIntegrandXi(real_t xi, void *par);
        struct intPParams {real_t pmax; real_t ntot; real_t mag; int qag_key; gsl_integration_workspace *gslws2; gsl_integration_workspace *gslws3;};
        static real_t analyticIntegrandP(real_t p, void *par);
    public:
        AvalancheSourceCH(const std::string& s) : UnitTest(s) {}

        DREAM::FVM::UnknownQuantityHandler *GetUnknownHandler(DREAM::FVM::Grid *g_fl, DREAM::FVM::Grid *g_hot, DREAM::FVM::Grid *g_re, const real_t n_re, const real_t n_tot, const real_t E_field);

        bool CheckConservativityCylindrical();
        bool CheckConservativityGeneralAnalytic();
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_DREAM_AVALANCHE_SOURCE_CH_HPP*/
