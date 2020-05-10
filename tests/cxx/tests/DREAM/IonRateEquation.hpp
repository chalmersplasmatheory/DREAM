#ifndef _DREAMTESTS_DREAM_ION_RATE_EQUATION_HPP
#define _DREAMTESTS_DREAM_ION_RATE_EQUATION_HPP

#include <string>
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
    class IonRateEquation : public UnitTest {
    private:
        len_t id_ions;

    public:
        IonRateEquation(const std::string& s) : UnitTest(s) {}

        DREAM::IonHandler *GetIonHandler(DREAM::FVM::Grid*, DREAM::FVM::UnknownQuantityHandler*);
        DREAM::FVM::UnknownQuantityHandler *GetUnknownHandler(DREAM::FVM::Grid*);

        bool CheckConservativity();
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_DREAM_ION_RATE_EQUATION_HPP*/
