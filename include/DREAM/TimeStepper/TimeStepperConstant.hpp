#ifndef _DREAM_TIME_STEPPER_CONSTANT_HPP
#define _DREAM_TIME_STEPPER_CONSTANT_HPP

#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class TimeStepperConstant : public TimeStepper {
    private:
        real_t dt;
        real_t tMax;

    public:
        TimeStepperConstant(const real_t tMax, const real_t dt, FVM::UnknownQuantityHandler *u)
            : TimeStepper(u), dt(dt), tMax(tMax) {}

        TimeStepperConstant(const real_t tMax, const len_t nt, FVM::UnknownQuantityHandler *u)
            : TimeStepper(u), tMax(tMax) {
            
            this->dt = tMax / nt;
        }

        virtual bool IsFinished(const real_t t) override { return (t >= tMax); }
        virtual real_t NextStep(const real_t) { return dt; }
    };
}

#endif/*_DREAM_TIME_STEPPER_CONSTANT_HPP*/
