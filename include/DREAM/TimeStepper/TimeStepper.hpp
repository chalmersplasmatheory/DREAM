#ifndef _DREAM_TIME_STEPPER_HPP
#define _DREAM_TIME_STEPPER_HPP

namespace DREAM { class EquationSystem; class TimeStepper; }

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class TimeStepper {
    public:
        typedef bool (*timestep_terminate_func_t)(void*, Simulation*);

    protected:
        FVM::UnknownQuantityHandler *unknowns;
        EquationSystem *eqsys;

        // Pointer to Solver object used for inverting equation system
        Solver *solver;

        // Python termination members
        timestep_terminate_func_t python_caller;
        void *python_terminate_func;

    public:
        TimeStepper(FVM::UnknownQuantityHandler *u, EquationSystem *eqsys)
            : unknowns(u), eqsys(eqsys) {}
        virtual ~TimeStepper() {}

        virtual real_t CurrentTime() const = 0;
        virtual void HandleException(FVM::FVMException&) = 0;
        virtual bool IsFinished() = 0;
        virtual bool IsSaveStep() = 0;
        virtual real_t MaxTime() const = 0;
        virtual real_t NextTime() = 0;
        virtual void PrintProgress() = 0;
        virtual void ValidateStep() = 0;

        bool PythonIsTerminate();

        void SetPythonCaller(timestep_terminate_func_t f) { this->python_caller = f; }
        void SetPythonTerminateFunc(void *f) { this->python_terminate_func = f; }
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
