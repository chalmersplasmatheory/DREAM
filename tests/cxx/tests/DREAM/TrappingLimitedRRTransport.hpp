#ifndef _DREAMTESTS_TRAPPING_LIMITED_RR_TRANSPORT_HPP
#define _DREAMTESTS_TRAPPING_LIMITED_RR_TRANSPORT_HPP

#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
    class TrappingLimitedRRTransport : public UnitTest {
    public:
        TrappingLimitedRRTransport(const std::string& name) : UnitTest(name) {}
        bool CheckGrid(DREAM::FVM::Grid*);
        virtual bool Run(bool) override;
    };
}
#endif
