/**
 * Implementation of the prescribed time stepper module.
 */

#include <cmath>

#include "DREAM/TimeStepper/TimeStepperPrescribed.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
TimeStepperPrescribed::TimeStepperPrescribed(
    const len_t ntimes, const real_t *t,
    FVM::UnknownQuantityHandler *u,
    EquationSystem *eqsys,
    const len_t nSaveSteps
) : TimeStepper(u, eqsys), nSaveSteps(nSaveSteps) {

    // Copy times to internal storage (ownership/lifetime safety).
    this->times.assign(t, t + ntimes);

    this->t0   = this->times.front();
    this->tMax = this->times.back();

    // Number of steps is one less than number of time points.
    this->Nt = (ntimes > 0 ? (ntimes - 1) : 0);

    InitSaveSteps();
}


/**
 * Returns the time of the most recently _completed_ time step.
 */
real_t TimeStepperPrescribed::CurrentTime() const {
    if (this->tIndex == 0)
        return this->t0;
    else
        return this->times[this->tIndex-1];
}


/**
 * Initialize the save-step checker.
 */
void TimeStepperPrescribed::InitSaveSteps() {
    if (this->nSaveSteps == 0) {
        return;
    } else if (this->nSaveSteps > this->Nt) {
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
bool TimeStepperPrescribed::IsFinished() {
    bool v = (this->tIndex >= this->Nt);
#ifdef DREAM_IS_PYTHON_LIBRARY
    return (v || this->PythonIsTerminate());
#else
    return v;
#endif
}


/**
 * Returns 'true' if the current time step should be saved to
 * the final output.
 */
bool TimeStepperPrescribed::IsSaveStep() {
    if (this->nSaveSteps == 0)
        return true;

    if (this->nextSaveStep_l == this->tIndex)
        return true;
    else
        return false;
}


/**
 * Returns the maximum time for this simulation.
 */
real_t TimeStepperPrescribed::MaxTime() const {
    return this->tMax;
}


/**
 * Returns the time of the next step.
 */
real_t TimeStepperPrescribed::NextTime() {
    this->tIndex++;

    // Update save step?
    if (this->nextSaveStep_l < this->tIndex) {
        this->nextSaveStep += this->dSaveStep;
        this->nextSaveStep_l = round(this->nextSaveStep);
    }

    return this->times[this->tIndex];
}


/**
 * Print current progress to stdout.
 */
void TimeStepperPrescribed::PrintProgress() {
    if (IsSaveStep())
        std::cout << "\x1B[1;32m" << this->tIndex << "\x1B[0m... ";
    else
        std::cout << this->tIndex << "... ";

    if (this->tIndex % 10 == 0) std::cout << std::endl;
}


/**
 * Validate the most recently taken time step.
 * (no validation needed for the prescribed time stepper...)
 */
void TimeStepperPrescribed::ValidateStep() {
}