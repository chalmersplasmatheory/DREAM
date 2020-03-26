#ifndef _DREAMTESTS_FVM_GRID_HPP
#define _DREAMTESTS_FVM_GRID_HPP

#include <string>
#include "FVM/Grid/Grid.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class Grid : public UnitTest {
    public:
        Grid(const std::string& name) : UnitTest(name) {}

        bool CheckGridRCylPXi();

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_GRID_HPP*/
