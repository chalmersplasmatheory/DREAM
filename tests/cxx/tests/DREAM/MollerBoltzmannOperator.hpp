#ifndef _DREAMTESTS_DREAM_BOLTZMANN_OPERATOR_HPP
#define _DREAMTESTS_DREAM_BOLTZMANN_OPERATOR_HPP

#include <string>

#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {

class MollerBoltzmannOperator : public UnitTest {
   public:
    MollerBoltzmannOperator(const std::string &s) : UnitTest(s) {}

    bool CheckBO_LinearityInFPrimary();
    bool CheckBO_JacobianFiniteDifferenceNt();
    bool CheckBO_GlobalProductionIdentity();
    bool CheckBO_HotRunawayGlobalProductionIdentity();
    bool CheckBO_TimeCachingRegression();
    bool CheckBO_NonNegativity();
    bool CheckBO_RadiusLocality();

    virtual bool Run(bool) override;
};

}  // namespace DREAMTESTS::_DREAM

#endif
