#ifndef _DREAM_TIME_STEPPER_CONSTANT_HPP
#define _DREAM_TIME_STEPPER_CONSTANT_HPP

#include <iostream>
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class TimeStepperConstant : public TimeStepper {
    private:
        len_t tIndex = 0;
        real_t dt;
        real_t t0 = 0;
        real_t tMax;
        len_t Nt;

    public:
        TimeStepperConstant(const real_t, const real_t, FVM::UnknownQuantityHandler*);
        TimeStepperConstant(const real_t, const len_t, FVM::UnknownQuantityHandler*);

        bool CheckNegative(const std::string&);

        virtual real_t CurrentTime() const override;
        virtual void HandleException(FVM::FVMException&) override;
        virtual bool IsFinished() override;
        virtual bool IsSaveStep() override;
        virtual real_t NextTime() override;
        virtual void PrintProgress() override;
        virtual void ValidateStep() override;
    };
}

#endif/*_DREAM_TIME_STEPPER_CONSTANT_HPP*/
