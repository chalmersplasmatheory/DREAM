#ifndef _DREAMTESTS_DREAM_KNOCK_ON_QUADRATURE_HPP
#define _DREAMTESTS_DREAM_KNOCK_ON_QUADRATURE_HPP
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
class KnockOnQuadrature : public UnitTest {
public:
    KnockOnQuadrature(const std::string& name) : UnitTest(name) {}
    bool Run(bool) override;
    bool CheckContactLevels();
    bool CheckPrimaryPitchSplits();
    bool CheckGeometry();
    bool CheckGaussRule();
    bool CheckPoloidalPanels();
    bool CheckValidation();
};
}
#endif
