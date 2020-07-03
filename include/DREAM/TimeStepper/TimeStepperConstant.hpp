#ifndef _DREAM_TIME_STEPPER_CONSTANT_HPP
#define _DREAM_TIME_STEPPER_CONSTANT_HPP

#include "DREAM/TimeStepper/TimeStepper.hpp"
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
        TimeStepperConstant(const real_t tMax, const real_t dt, FVM::UnknownQuantityHandler *u)
            : TimeStepper(u), dt(dt), tMax(tMax) { this->Nt = round(tMax/dt); }

        TimeStepperConstant(const real_t tMax, const len_t nt, FVM::UnknownQuantityHandler *u)
            : TimeStepper(u), tMax(tMax), Nt(nt) {
            
            this->dt = tMax / nt;
        }

        virtual real_t CurrentTime() const override { return (this->t0 + this->tIndex*this->dt); }
        virtual bool IsFinished() override { return (this->tIndex>=this->Nt); }
        virtual real_t NextTime() {
            this->tIndex++;
            return CurrentTime();
        }
    };
}

#endif/*_DREAM_TIME_STEPPER_CONSTANT_HPP*/
