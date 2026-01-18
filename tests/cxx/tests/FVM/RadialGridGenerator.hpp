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
    bool CheckNumericBGeometryCalculations();
    bool CheckNumericBAgreesWithAnalyticB();

    virtual bool Run(bool) override;

   private:
    bool TestGeometricQuantitiesAtTheta(
        len_t ir, real_t theta, DREAM::FVM::RadialGrid* rg, real_t atol, real_t rtol
    );
    void printComparison(
        real_t numeric, real_t analytic, const std::string& name1, const std::string& name2,
        real_t atol, real_t rtol, len_t n_indent = 4
    );
};
}  // namespace DREAMTESTS::FVM

#endif /*_DREAMTESTS_FVM_RADIAL_GRID_GENERATOR_HPP*/
