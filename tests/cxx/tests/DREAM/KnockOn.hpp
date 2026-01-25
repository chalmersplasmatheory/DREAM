#ifndef _DREAMTESTS_DREAM_KNOCK_ON_HPP
#define _DREAMTESTS_DREAM_KNOCK_ON_HPP

#include <string>

#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
class KnockOn : public UnitTest {
   public:
    KnockOn(const std::string &s) : UnitTest(s) {}

    bool CheckDeltaMirrorProperties();
    bool CheckDeltaQuadratureConvergence();
    bool CheckLocalContributionConservationProperty();
    bool CheckDeltaConservationProperty();
    bool CheckDeltaConservationPropertyDifferentGrids();
    bool CheckAgreementWithOldRPTerm();
    bool CheckCylindricalDeltaCalculation();
    bool CheckMollerFluxIntegration();
    bool CheckMollerDifferentialConvergesToInfiniteLimit();
    bool CheckMollerFluxConvergesToInfiniteLimit();
    bool CheckXiStarConvergesToInfiniteLimit();

    virtual bool Run(bool) override;

   private:
    bool _checkLocalContributionConservationProperty(
        const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary, real_t tol
    );
    bool _checkDeltaConservationProperty(
        const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary, real_t tol
    );
};
}  // namespace DREAMTESTS::_DREAM

#endif /*_DREAMTESTS_DREAM_KNOCK_ON_HPP*/
