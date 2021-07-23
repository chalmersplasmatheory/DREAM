/**
 * Access to details from the 'Simulation' object.
 */

#include "DREAM/Simulation.hpp"
#include "pyface/simulation.hpp"


DREAM::Simulation *get_simulation_from_capsule(PyObject *args) {
    PyObject *sim_capsule;
    if (!PyArg_ParseTuple(args, "O", &sim_capsule) || !PyCapsule_CheckExact(sim_capsule)) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "Expected argument to be a Python capsule."
        );
        return NULL;
    }

    return reinterpret_cast<DREAM::Simulation*>(
        PyCapsule_GetPointer(sim_capsule, "sim")
    );
}


extern "C" {

/**
 * Returns the current time from the given C++ Simulation object.
 *
 * PYTHON PARAMETERS
 * sim_capsule: PyCapsule object containing a pointer to the
 *              DREAM::Simulation object to access.
 */
static PyObject *dreampy_get_current_time(
    PyObject* /*self*/, PyObject *args
) {
    DREAM::Simulation *sim = get_simulation_from_capsule(args);

    if (sim == NULL)
        return NULL;

    return PyFloat_FromDouble(sim->GetEquationSystem()->GetCurrentTime());
}

/**
 * Returns the maximum simulation time from the given C++ Simulation object.
 *
 * PYTHON PARAMETERS
 * sim_capsule: PyCapsule object containing a pointer to the
 *              DREAM::Simulation object to access.
 */
static PyObject *dreampy_get_max_time(
    PyObject* /*self*/, PyObject *args
) {
    DREAM::Simulation *sim = get_simulation_from_capsule(args);

    if (sim == NULL)
        return NULL;

    return PyFloat_FromDouble(sim->GetEquationSystem()->GetMaxTime());
}

}

