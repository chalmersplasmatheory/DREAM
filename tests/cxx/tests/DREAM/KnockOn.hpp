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
    bool CheckPoloidalNumerators(bool extended=false);
    bool CheckPitchCellAverages(bool extended=false);
    bool CheckLocalContributionConservationProperty();
    bool CheckDeltaConservationProperty();
    bool CheckDeltaConservationPropertyConstantB();
    bool CheckDeltaConservationPropertyDifferentGrids();
    bool CheckAgreementWithOldRPTerm();
    bool CheckCylindricalDeltaCalculation();
    bool CheckMollerFluxIntegration();
    bool CheckMollerDifferentialConvergesToInfiniteLimit();
    bool CheckMollerFluxConvergesToInfiniteLimit();
    bool CheckXiStarConvergesToInfiniteLimit();
    bool CheckMollerSConservationProperty();

    virtual bool Run(bool) override;
    bool RunExtended();

   private:
    bool _checkLocalContributionConservationProperty(
        const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary, real_t tol
    );
    bool _checkDeltaConservationProperty(
        const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary, real_t tol
    );
    bool _checkMollerSConservationProperty(
        const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary
    );
};
}  // namespace DREAMTESTS::_DREAM

#endif /*_DREAMTESTS_DREAM_KNOCK_ON_HPP*/
