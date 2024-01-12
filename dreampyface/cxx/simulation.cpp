/**
 * Access to details from the 'Simulation' object.
 */

#include <vector>
#include "DREAM/Simulation.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "pyface/numpy.h"
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

/**
 * Returns the radial grid used for the simulation.
 *
 * PYTHON PARAMETERS
 * sim: PyCapsule object containing a pointer to the DREAM::Simulation
 *      object to access.
 */
static PyObject *dreampy_get_radius_vector(
    PyObject* /*self*/, PyObject *args
) {
    DREAM::Simulation *sim = get_simulation_from_capsule(args);

    if (sim == NULL)
        return NULL;

    DREAM::FVM::RadialGrid *rg = sim->GetEquationSystem()->GetFluidGrid()->GetRadialGrid();
    const real_t *radii = rg->GetR();

    npy_intp l = rg->GetNr();
    PyObject *arr = PyArray_SimpleNew(1, &l, NPY_DOUBLE);
    real_t *p = reinterpret_cast<real_t*>(
        PyArray_DATA(reinterpret_cast<PyArrayObject*>(arr))
    );

    for (npy_intp i = 0; i < l; i++)
        p[i] = radii[i];

    return arr;
}

/**
 * Returns the list of time steps taken so far.
 *
 * PYTHON PARAMETERS
 * sim: PyCapsule object containing a pointer to the DREAM::Simulation
 *      object to access.
 */
static PyObject *dreampy_get_time_vector(
    PyObject* /*self*/, PyObject *args
) {
    DREAM::Simulation *sim = get_simulation_from_capsule(args);

    if (sim == NULL)
        return NULL;

    std::vector<real_t>& times = sim->GetEquationSystem()->GetTimes();

    npy_intp l = times.size();
    PyObject *arr = PyArray_SimpleNew(1, &l, NPY_DOUBLE);
    real_t *p = reinterpret_cast<real_t*>(
        PyArray_DATA(reinterpret_cast<PyArrayObject*>(arr))
    );

    for (npy_intp i = 0; i < l; i++)
        p[i] = times[i];

    return arr;

}

}

