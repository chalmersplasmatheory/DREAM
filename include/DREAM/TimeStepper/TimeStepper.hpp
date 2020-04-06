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

        virtual bool IsFinished(const real_t currentTime) = 0;
        virtual real_t NextStep(const real_t currenTime) = 0;
    };
}

#endif/*_DREAM_TIME_STEPPER_HPP*/
