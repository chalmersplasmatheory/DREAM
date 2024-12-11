#ifndef _DREAMTESTS_DREAM_AVALANCHE_SOURCE_CH_HPP
#define _DREAMTESTS_DREAM_AVALANCHE_SOURCE_CH_HPP

#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
    class AvalancheSourceCH : public UnitTest {
    public:
        AvalancheSourceCH(const std::string& s) : UnitTest(s) {}

        DREAM::FVM::UnknownQuantityHandler *GetUnknownHandler(DREAM::FVM::Grid*,const real_t f_re,const real_t f_hot, const real_t n_tot, const real_t E_field);

        bool CheckConservativityCylindrical();
        //bool CheckConservativityGeneralAnalytic();
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_DREAM_AVALANCHE_SOURCE_CH_HPP*/
