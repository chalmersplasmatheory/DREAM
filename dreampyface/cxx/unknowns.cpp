/**
 * Functions for obtaining information and data for
 * unknown quantities in the simulation.
 */

#include "DREAM/Simulation.hpp"
#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "pyface/unknowns.hpp"

extern "C" {

/**
 * Get a dictionary listing all the unknowns of the equation
 * system, their size, as well as a description of their definition
 * (roughly equivalent to the equation system description printed
 * at the start of a DREAM simulation)
 */
static PyObject *dreampy_get_unknowns(
    PyObject* /*self*/, PyObject *args
) {
    DREAM::Simulation *sim = get_simulation_from_capsule(args);

    if (sim == NULL)
        return NULL;

    PyObject *dict = PyDict_New();
    for (auto u : sim->GetEquationSystem()->GetUnknownHandler()->GetUnknowns()) {
        PyObject *inner = PyDict_New();

        PyDict_SetItemString(inner, "description", PyUnicode_FromString(u->GetDescription().c_str()));
        PyDict_SetItemString(inner, "equation", PyUnicode_FromString(u->GetEquationDescription().c_str()));
        PyDict_SetItemString(inner, "nelements", PyLong_FromUnsignedLongLong(u->NumberOfElements()));
        PyDict_SetItemString(inner, "nmultiples", PyLong_FromUnsignedLongLong(u->NumberOfMultiples()));

        PyDict_SetItemString(dict, u->GetName().c_str(), inner);
    }

    return dict;
}

/**
 * Return info about a single named unknown quantity.
 *
 * PYTHON PARAMETERS
 * sim:  Pointer to C++ Simulation object.
 * name: Name of unknown quantity to obtain information for.
 */
static PyObject *dreampy_get_unknown_info(
    PyObject* /*self*/, PyObject *args, PyObject *kwargs
) {
    static const char *kwlist[] = {"simulation", "name", NULL};
    PyObject *simulation;
    const char *name;

    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, "Os", const_cast<char**>(kwlist),
        &simulation, &name
    )) {
        return NULL;
    }

    DREAM::Simulation *sim = reinterpret_cast<DREAM::Simulation*>(
        PyCapsule_GetPointer(simulation, "sim")
    );

    DREAM::FVM::UnknownQuantityHandler *uqn = sim->GetEquationSystem()->GetUnknownHandler();
    len_t id;
    bool exists = true;
    try {
        id = uqn->GetUnknownID(name);
    } catch (DREAM::FVM::FVMException& ex) {
        // Convert to Python exception
        PyErr_SetString(
            PyExc_RuntimeError,
            ex.what()
        );
        exists = false;
    }

    if (!exists)
        return NULL;

    DREAM::FVM::UnknownQuantity *u = uqn->GetUnknown(id);

    PyObject *dct = PyDict_New();

    PyDict_SetItemString(dct, "description", PyUnicode_FromString(u->GetDescription().c_str()));
    PyDict_SetItemString(dct, "equation", PyUnicode_FromString(u->GetEquationDescription().c_str()));
    PyDict_SetItemString(dct, "nelements", PyLong_FromUnsignedLongLong(u->NumberOfElements()));
    PyDict_SetItemString(dct, "nmultiples", PyLong_FromUnsignedLongLong(u->NumberOfMultiples()));

    return dct;
}

}

