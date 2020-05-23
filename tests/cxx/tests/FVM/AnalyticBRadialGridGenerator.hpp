#ifndef _DREAMTESTS_FVM_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP
#define _DREAMTESTS_FVM_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP

#include "UnitTest.hpp"
//#include "FVM/Grid/RadialGrid.hpp"

namespace DREAMTESTS::FVM {
    class AnalyticBRadialGridGenerator : public UnitTest {
    public:
        AnalyticBRadialGridGenerator(const std::string& name) : UnitTest(name){}

        DREAM::FVM::Grid *grid;

        virtual void Initialize();

        virtual bool TestGeneralBounceAverage();
        virtual bool TestGeneralFluxSurfaceAverage();
        virtual bool CompareBounceAverageMethods();
        virtual bool Run(bool) override;

    };
}

#endif/*_DREAMTESTS_FVM_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP*/
