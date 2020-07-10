#ifndef _DREAMTESTS_FVM_PXI_EXTERNAL_KINETIC_KINETIC_HPP
#define _DREAMTESTS_FVM_PXI_EXTERNAL_KINETIC_KINETIC_HPP

#include <string>
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::FVM {
    class PXiExternalKineticKinetic : public UnitTest {
    private:
    public:
        PXiExternalKineticKinetic(const std::string& s) : UnitTest(s) {}

        bool CheckConsistency();
        bool CompareToPXiExternalLoss();

        bool Check(
            bool (PXiExternalKineticKinetic::*)(
                DREAM::FVM::Operator*, const std::string&,
                DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*
            ),
            len_t nxi_re=0, bool sameSizeRE=false
        );
        bool CheckPXiExternalLoss(
            DREAM::FVM::Operator*, const std::string&,
            DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*
        );
        bool CheckConservativity(
            DREAM::FVM::Operator*, const std::string&,
            DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*
        );

        real_t *ConvertFlux(const real_t*, DREAM::FVM::Grid*, DREAM::FVM::Grid*);

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_FVM_PXI_EXTERNAL_KINETIC_KINETIC_HPP*/
