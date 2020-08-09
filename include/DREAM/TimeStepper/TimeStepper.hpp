#ifndef _DREAM_TIME_STEPPER_HPP
#define _DREAM_TIME_STEPPER_HPP

#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class TimeStepper {
    protected:
        FVM::UnknownQuantityHandler *unknowns;

    public:
        TimeStepper(FVM::UnknownQuantityHandler *u)
            : unknowns(u) {}
        virtual ~TimeStepper() {}

        virtual real_t CurrentTime() const = 0;
        virtual bool IsFinished() = 0;
        virtual real_t NextTime() = 0;
        virtual void PrintProgress() = 0;
        virtual bool IsSaveStep() = 0;
    };

    class TimeStepperException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        TimeStepperException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("TimeStepper");
        }
    };
}

#endif/*_DREAM_TIME_STEPPER_HPP*/
