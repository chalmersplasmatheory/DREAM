/**
 * Implementation of the constant time stepper module.
 */

#include "DREAM/IO.hpp"
#include "DREAM/TimeStepper/TimeStepperConstant.hpp"


using namespace DREAM;


/**
 * Constructors.
 */
TimeStepperConstant::TimeStepperConstant(
    const real_t tMax, const real_t dt, FVM::UnknownQuantityHandler *u,
    const len_t nSaveSteps
) : TimeStepper(u), dt(dt), tMax(tMax), nSaveSteps(nSaveSteps) {

    this->Nt = round(tMax/dt);
    InitSaveSteps();
}
TimeStepperConstant::TimeStepperConstant(
    const real_t tMax, const len_t nt, FVM::UnknownQuantityHandler *u,
    const len_t nSaveSteps
) : TimeStepper(u), tMax(tMax), Nt(nt), nSaveSteps(nSaveSteps) {
    
    this->dt = tMax / nt;
    InitSaveSteps();
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
 * Initialize the save-step checker.
 */
void TimeStepperConstant::InitSaveSteps() {
    if (this->nSaveSteps == 0) return;
    else if (this->nSaveSteps > this->Nt) {
        this->nSaveSteps = 0;
        return;
    }

    this->dSaveStep = real_t(this->Nt) / real_t(this->nSaveSteps);
    this->nextSaveStep = this->dSaveStep;
    this->nextSaveStep_l = round(this->nextSaveStep);   // >= 1
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
    if (this->nSaveSteps == 0)
        return true;
    
    if (this->nextSaveStep_l == this->tIndex)
        return true;
    else
        return false;
}

/**
 * Returns the time of the next step.
 */
real_t TimeStepperConstant::NextTime() {
    this->tIndex++;

    // Update save step?
    if (this->nextSaveStep_l < this->tIndex) {
        this->nextSaveStep += this->dSaveStep;
        this->nextSaveStep_l = round(this->nextSaveStep);
    }

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
