#ifndef _DREAMTESTS_DREAM_MOLLER_DELTA_ANGLE_KERNEL_HPP
#define _DREAMTESTS_DREAM_MOLLER_DELTA_ANGLE_KERNEL_HPP

#include <string>

#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {

class MollerDeltaAngleKernel : public UnitTest {
   public:
    MollerDeltaAngleKernel(const std::string &s) : UnitTest(s) {}

    bool CheckMDK_DeltaInterpolationConservation();
    bool CheckMDK_KinematicEdgeCases();
    bool CheckMDK_IsotropicAndDipoleMoments();
    bool CheckMDK_GridRebuilt();
    virtual bool Run(bool) override;
};

}  // namespace DREAMTESTS::_DREAM

#endif
