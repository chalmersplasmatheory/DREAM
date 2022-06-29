#ifndef _DREAM_TIME_STEPPER_HPP
#define _DREAM_TIME_STEPPER_HPP

#include "DREAM/Solver/Solver.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class TimeStepper {
    protected:
        FVM::UnknownQuantityHandler *unknowns;

        // Pointer to Solver object used for inverting equation system
        Solver *solver;

    public:
        TimeStepper(FVM::UnknownQuantityHandler *u)
            : unknowns(u) {}
        virtual ~TimeStepper() {}

        virtual bool CheckNegative(const std::string&);
        virtual void HandleException(FVM::FVMException&);

        virtual real_t CurrentTime() const = 0;
        virtual bool IsFinished() = 0;
        virtual bool IsSaveStep() = 0;
        virtual real_t NextTime() = 0;
        virtual void PrintProgress() = 0;
        virtual void ValidateStep() = 0;

        void SetSolver(Solver *s) { this->solver = s; }
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
