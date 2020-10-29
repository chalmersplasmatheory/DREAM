/**
 * Implementation of the constant time stepper module.
 */

#include <vector>
#include "DREAM/IO.hpp"
#include "DREAM/TimeStepper/TimeStepperConstant.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructors.
 */
TimeStepperConstant::TimeStepperConstant(
    const real_t tMax, const real_t dt, FVM::UnknownQuantityHandler *u,
    vector<len_t>& nontrivials, bool useBackupInverter
) : TimeStepper(u, nontrivials), dt(dt), tMax(tMax), useBackupInverter(useBackupInverter) {
    this->Nt = round(tMax/dt);
}
TimeStepperConstant::TimeStepperConstant(
    const real_t tMax, const len_t nt, FVM::UnknownQuantityHandler *u,
    vector<len_t>& nontrivials, bool useBackupInverter
) : TimeStepper(u, nontrivials), tMax(tMax), Nt(nt), useBackupInverter(useBackupInverter) {
    this->dt = tMax / nt;
}

/**
 * Returns the time of the most recently _completed_ time step.
 */
real_t TimeStepperConstant::CurrentTime() const {
    if (this->tIndex == 0)
        return this->t0;
    else
        return (this->t0 + (this->tIndex-1)*this->dt);
}

/**
 * This method is called when an exception was thrown while
 * the next time step was being taken. In the constant time
 * stepper, we choose to forward this exception and just indicate
 * to the user that it might be due to the time step being too
 * short.
 *
 * ex: The exception that was caught.
 */
void TimeStepperConstant::HandleException(FVM::FVMException &ex) {
    // It has been observed that some linear solvers sometimes fail where
    // others succeed. Thus we can try to complete the step by temporarily
    // switching to the "backup" solver.
    if (!this->useBackupInverter || solver->IsUsingBackupInverter() || !solver->HasBackupInverter()) {
        DREAM::IO::PrintError("TimeStepper: Exception caught during time stepping.");
        DREAM::IO::PrintError("TimeStepper: Perhaps the exception could be avoided by decreasing the time step length?");
        throw ex;
    } else {
        DREAM::IO::PrintInfo("Switching to backup inverter");
        solver->SetUseBackupInverter(true);

        RestoreInitialSolution(1, false);
    }
}

/**
 * Returns 'true' if the time stepper has reached the maximum time.
 */
bool TimeStepperConstant::IsFinished() {
    return (this->tIndex>=this->Nt);
}

/**
 * Returns 'true' if the current time step should be saved to
 * the final output. (Currently, we save all time steps)
 */
bool TimeStepperConstant::IsSaveStep() {
    return true;
}

/**
 * Returns the time of the next step.
 */
real_t TimeStepperConstant::NextTime() {
    this->tIndex++;
    this->initTime = this->t0;
    return (this->t0 + this->tIndex*this->dt);
}

/**
 * Print current progress to stdout.
 */
void TimeStepperConstant::PrintProgress() {
    if (IsSaveStep())
        std::cout << "\x1B[1;32m" << this->tIndex << "\x1B[0m... ";
    else
        std::cout << this->tIndex << "... ";

    if (this->tIndex % 10 == 0) std::cout << std::endl;
}

/**
 * Validate the most recently taken time step.
 * (no validation needed for the constant time stepper...)
 */
void TimeStepperConstant::ValidateStep() {
}
