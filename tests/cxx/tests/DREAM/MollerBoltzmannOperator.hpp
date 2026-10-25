#ifndef _DREAMTESTS_DREAM_BOLTZMANN_OPERATOR_HPP
#define _DREAMTESTS_DREAM_BOLTZMANN_OPERATOR_HPP

#include <string>
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
class MollerBoltzmannOperator : public UnitTest {
   public:
    MollerBoltzmannOperator(const std::string &s) : UnitTest(s) {}
    virtual bool Run(bool) override;

   private:
    struct Fixture;
    bool CheckBO_LinearityInFPrimary(const Fixture&);
    bool CheckBO_JacobianFiniteDifferenceNt(const Fixture&);
    bool CheckBO_GlobalProductionIdentity(const Fixture&);
    bool CheckBO_HotRunawayGlobalProductionIdentity();
    bool CheckBO_TimeCachingRegression(const Fixture&);
    bool CheckBO_MatrixAssemblyAndFluidIdentity(const Fixture&);
    bool CheckBO_NonNegativity(const Fixture&);
    bool CheckBO_RadiusLocality(const Fixture&);
    bool CheckBO_IsotropicZeroCurrentAndReflection();
};
} // namespace DREAMTESTS::_DREAM
#endif
