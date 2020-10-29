#ifndef _DREAM_TIME_STEPPER_HPP
#define _DREAM_TIME_STEPPER_HPP

#include "DREAM/Solver/Solver.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class TimeStepper {
    protected:
        FVM::UnknownQuantityHandler *unknowns;

        // List of non-trivial unknowns
        std::vector<len_t> nontrivials;

        // Pointer to Solver object used for inverting equation system
        Solver *solver;

        real_t initTime = 0;
        len_t sol_size = 0;
        // Initial solution (before taking the first half step)
        real_t *sol_init=nullptr;

    public:
        TimeStepper(FVM::UnknownQuantityHandler *u, std::vector<len_t>& nt)
            : unknowns(u), nontrivials(nt) {}
        virtual ~TimeStepper() {}

        virtual void AllocateSolutions(const len_t);
        virtual void DeallocateSolutions();
        void RestoreInitialSolution(const len_t, bool pushinit=true);

        virtual real_t CurrentTime() const = 0;
        virtual void HandleException(FVM::FVMException&) = 0;
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
