#ifndef _DREAMTESTS_EXTERNAL_BOUNDARY_HPP
#define _DREAMTESTS_EXTERNAL_BOUNDARY_HPP

#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
class ExternalBoundary : public UnitTest {
public:
    ExternalBoundary(const std::string& name) : UnitTest(name) {}
    bool Run(bool) override;
private:
    bool CheckGrid(DREAM::FVM::Grid *, DREAM::FVM::Grid *);
};
}
#endif
