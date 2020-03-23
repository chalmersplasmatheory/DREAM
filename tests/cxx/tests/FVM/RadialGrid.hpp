#ifndef _DREAMTESTS_FVM_RADIAL_GRID_HPP
#define _DREAMTESTS_FVM_RADIAL_GRID_HPP

#include <string>
#include "FVM/Grid/RadialGrid.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class RadialGrid : public UnitTest {
    public:
        RadialGrid(const std::string& name) : UnitTest(name) {}

        bool CheckGeneralGrid();

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_RADIAL_GRID_HPP*/
