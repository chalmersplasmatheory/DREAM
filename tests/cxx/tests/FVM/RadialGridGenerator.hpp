#ifndef _DREAMTESTS_FVM_RADIAL_GRID_GENERATOR_HPP
#define _DREAMTESTS_FVM_RADIAL_GRID_GENERATOR_HPP

#include <string>
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class RadialGridGenerator : public UnitTest {

    public:
        RadialGridGenerator(const std::string& name) : UnitTest(name) {}

        bool CheckAnalyticBGeometryCalculations();

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_RADIAL_GRID_GENERATOR_HPP*/
