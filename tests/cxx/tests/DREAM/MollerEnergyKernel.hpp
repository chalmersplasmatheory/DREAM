#ifndef _DREAMTESTS_DREAM_MOLLER_ENERGY_KERNEL_HPP
#define _DREAMTESTS_DREAM_MOLLER_ENERGY_KERNEL_HPP

#include <string>

#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {

class MollerEnergyKernel : public UnitTest {
   public:
    MollerEnergyKernel(const std::string &s) : UnitTest(s) {}

    bool CheckMEK_MollerKernelConservation();

    virtual bool Run(bool) override;
};

}  // namespace DREAMTESTS::_DREAM

#endif
