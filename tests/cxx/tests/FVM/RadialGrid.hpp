#ifndef _TQSTESTS_FVM_RADIAL_GRID_HPP
#define _TQSTESTS_FVM_RADIAL_GRID_HPP

#include <string>
#include "FVM/Grid/RadialGrid.hpp"
#include "UnitTest.hpp"

namespace TQSTESTS::FVM {
    class RadialGrid : public UnitTest {
    public:
        RadialGrid(const std::string& name) : UnitTest(name) {}

        bool CheckGeneralGrid();

        virtual bool Run(bool) override;
    };
}

#endif/*_TQSTESTS_FVM_RADIAL_GRID_HPP*/
