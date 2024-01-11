/**
 * Implementation of routines common for all TimeStepper classes.
 */

#include <string>
#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/TimeStepper/TimeStepper.hpp"


using namespace DREAM;


/**
 * Calls a previously provided Python routine to determine
 * whether or not time stepping should continue.
 */
bool TimeStepper::PythonIsTerminate() {
    if (this->python_terminate_func != nullptr && this->python_caller != nullptr)
        return this->python_caller(this->python_terminate_func, this->eqsys->GetSimulation());
    else
        return false;
}

/**
 * Check if a given unknown contains negative elements.
 */
bool TimeStepper::CheckNegative(const std::string& name) {
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
 * This method is called when an exception was thrown while
 * the next time step was being taken. In the constant time
 * stepper, we choose to forward this exception and just indicate
 * to the user that it might be due to the time step being too
 * short.
 *
 * ex: The exception that was caught.
 */
void TimeStepper::HandleException(FVM::FVMException &ex) {
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

