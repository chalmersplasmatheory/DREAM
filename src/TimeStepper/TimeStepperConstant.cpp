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
    EquationSystem *eqsys, const len_t nSaveSteps
) : TimeStepper(u, eqsys), dt(dt), tMax(tMax), nSaveSteps(nSaveSteps) {

    this->Nt = round(tMax/dt);
    InitSaveSteps();
}
TimeStepperConstant::TimeStepperConstant(
    const real_t tMax, const len_t nt, FVM::UnknownQuantityHandler *u,
    EquationSystem *eqsys, const len_t nSaveSteps
) : TimeStepper(u, eqsys), tMax(tMax), Nt(nt), nSaveSteps(nSaveSteps) {
    
    this->dt = tMax / nt;
    InitSaveSteps();
}

/**
 * Check if a given unknown contains negative elements.
 */
bool TimeStepperConstant::CheckNegative(const std::string& name) {
    len_t uqtyid = unknowns->GetUnknownID(name);

    if (unknowns->HasUnknown(name)) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(uqtyid);
        const real_t *data = uqty->GetData();
        for (len_t i = 0; i < uqty->NumberOfElements(); i++)
            if (data[i] < 0)
                return true;
    }

    return false;
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
    DREAM::IO::PrintError("TimeStepper: Exception caught during time stepping.");
    DREAM::IO::PrintError("TimeStepper: Perhaps the exception could be avoided by decreasing the time step length?");

//    if (CheckNegative(OptionConstants::UQTY_ION_SPECIES))
//        DREAM::IO::PrintError("TimeStepper: Ion density 'n_i' is negative.");
    if (CheckNegative(OptionConstants::UQTY_N_COLD))
        DREAM::IO::PrintError("TimeStepper: Cold electron density 'n_cold' is negative.");
    if (CheckNegative(OptionConstants::UQTY_T_COLD))
        DREAM::IO::PrintError("TimeStepper: Cold electron temperature 'T_cold' is negative.");

    throw ex;
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
    bool v = (this->tIndex>=this->Nt);
#ifdef DREAM_IS_PYTHON_LIBRARY
    return (v || this->PythonIsTerminate());
#endif
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
 * Returns the maximum time for this simulation.
 */
real_t TimeStepperConstant::MaxTime() const {
    return this->tMax;
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
