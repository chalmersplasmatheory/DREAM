#ifndef _DREAMTESTS_DREAM_KNOCK_ON_HPP
#define _DREAMTESTS_DREAM_KNOCK_ON_HPP

#include "UnitTest.hpp"
#include <string>

namespace DREAMTESTS::_DREAM {
class KnockOn : public UnitTest {
  public:
    KnockOn(const std::string &s) : UnitTest(s) {}

    bool CheckDeltaMirrorProperties();
    bool CheckDeltaQuadratureConvergence();
    bool CheckDeltaConservationProperty();
    bool CheckAgreementWithOldRPTerm();
    bool CheckCylindricalDeltaCalculation();

    virtual bool Run(bool) override;
};
} // namespace DREAMTESTS::_DREAM

#endif /*_DREAMTESTS_DREAM_KNOCK_ON_HPP*/
