/**
 * Implementation of routines common for all TimeStepper classes.
 */

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

