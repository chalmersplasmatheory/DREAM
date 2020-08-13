#ifndef _DREAM_TIME_STEPPER_CONSTANT_HPP
#define _DREAM_TIME_STEPPER_CONSTANT_HPP

#include <iostream>
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

        virtual real_t CurrentTime() const override {
            if (this->tIndex == 0) return this->t0;
            else
                return (this->t0 + (this->tIndex-1)*this->dt);
        }
        virtual bool IsFinished() override { return (this->tIndex>=this->Nt); }
        // Save all time steps
        virtual bool IsSaveStep() override { return true; }
        virtual real_t NextTime() override {
            this->tIndex++;
            return (this->t0 + this->tIndex*this->dt);
        }

        // Print current progress to stdout...
        virtual void PrintProgress() override {
            if (IsSaveStep())
                std::cout << "\x1B[1;32m" << this->tIndex << "\x1B[0m... ";
            else
                std::cout << this->tIndex << "... ";

            if (this->tIndex % 10 == 0) std::cout << std::endl;
        }

        // No validation needed for the constant time stepper...
        virtual void ValidateStep() override {}
    };
}

#endif/*_DREAM_TIME_STEPPER_CONSTANT_HPP*/
