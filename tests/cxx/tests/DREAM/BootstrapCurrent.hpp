#ifndef _DREAMTESTS_DREAM_BOOTSRTAP_CURRENT_HPP
#define _DREAMTESTS_DREAM_BOOTSTRAP_CURRENT_HPP

#include <string>
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
    class BootstrapCurrent : public UnitTest {
    private:

    public:
        BootstrapCurrent(const std::string& s) : UnitTest(s) {}

        DREAM::IonHandler *GetIonHandler(DREAM::FVM::Grid*, DREAM::FVM::UnknownQuantityHandler*);
        DREAM::FVM::UnknownQuantityHandler *GetUnknownHandler(DREAM::FVM::Grid*);

        bool CheckBootstrap();
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_DREAM_ION_RATE_EQUATION_HPP*/
