#ifndef _DREAMTESTS_DREAM_BOUNDARY_FLUX_HPP
#define _DREAMTESTS_DREAM_BOUNDARY_FLUX_HPP

#include <string>
#include "FVM/Equation/Equation.hpp"
#include "FVM/Grid/Grid.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
    class BoundaryFlux : public UnitTest {
    private:
    public:
        BoundaryFlux(const std::string& s) : UnitTest(s) {}

        bool CheckPXiToFluid();
        bool CheckPXiToFluid(
            const DREAM::FVM::Equation*, const std::string&,
            DREAM::FVM::Grid*, DREAM::FVM::Grid*
        );
        virtual bool Run(bool) override;

        bool VerifyVectorsAreOpposite(
            const std::string&, const std::string&,
            const len_t, const len_t,
            const real_t*, const real_t*
        );
    };
}

#endif/*_DREAMTESTS_DREAM_BOUNDARY_FLUX_HPP*/
