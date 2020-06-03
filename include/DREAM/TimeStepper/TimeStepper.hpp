#ifndef _DREAM_TIME_STEPPER_HPP
#define _DREAM_TIME_STEPPER_HPP

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
    };
}

#endif/*_DREAM_TIME_STEPPER_HPP*/
