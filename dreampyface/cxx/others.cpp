/**
 * Functions for obtaining information and data for
 * other quantities in the simulation.
 */

#include "DREAM/OtherQuantity.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/Simulation.hpp"
#include "FVM/FVMException.hpp"


/**
 * Returns the OtherQuantity with the specified name from the
 * given Simulation object. If 'NULL' is returned, this indicates
 * that a Python exception has been thrown and should be forwarded.
 *
 * sim:  Simulation object to fetch other quantity from.
 * name: Name of other quantity to obtain.
 */
DREAM::OtherQuantity *get_other(DREAM::Simulation *sim, const char *name) {
    DREAM::OtherQuantityHandler *oqn = sim->GetEquationSystem()->GetOtherQuantityHandler();

    DREAM::OtherQuantity *o = NULL;
    try {
        o = oqn->GetByName(name);

        if (o == NULL)
            PyErr_SetString(
                PyExc_RuntimeError,
                "Other quantity does not exist."
            );
    } catch (DREAM::FVM::FVMException& ex) {
        // Convert to Python exception
        PyErr_SetString(
            PyExc_RuntimeError,
            ex.what()
        );
    }

    return o;
}


extern "C" {

/**
 * Get a dictionary listing all the other quantities of the equation
 * system, their size, as well as a description of their definition.
 */
static PyObject *dreampy_get_others(
    PyObject* /*self*/, PyObject *args
) {
    DREAM::Simulation *sim = get_simulation_from_capsule(args);

    if (sim == NULL)
        return NULL;

    PyObject *dict = PyDict_New();
    for (auto u : sim->GetEquationSystem()->GetOtherQuantityHandler()->GetRegisteredQuantities()) {
        PyObject *inner = PyDict_New();

        PyDict_SetItemString(inner, "description", PyUnicode_FromString(u->GetDescription().c_str()));
        PyDict_SetItemString(inner, "nelements", PyLong_FromUnsignedLongLong(u->NumberOfElements()));
        PyDict_SetItemString(inner, "nmultiples", PyLong_FromUnsignedLongLong(u->NumberOfMultiples()));

        PyDict_SetItemString(dict, u->GetName().c_str(), inner);
    }

    return dict;
}

/**
 * Return the currently calculated data for the named other quantity.
 *
 * PYTHON PARAMETERS
 * sim:  Pointer to C++ Simulation object.
 * name: Name of unknown quantity to obtain information for.
 */
static PyObject *dreampy_get_other_data(
    PyObject* /*self*/, PyObject *args, PyObject *kwargs
) {
    static const char *kwlist[] = {"ptr", "name", NULL};
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

    DREAM::OtherQuantity *o = get_other(sim, name);

    if (o == NULL)
        return NULL;

    SFile_Python *sfp = new SFile_Python();
    o->SaveSFile(sfp);

    PyObject *dict = sfp->GetPythonDict();
    delete sfp;

    return dict;
}

/**
 * Return infor about a single named other quantity.
 *
 * PYTHON PARAMETERS
 * sim:  Pointer to C++ Simulation object.
 * name: Name of unknown quantity to obtain information for.
 */
static PyObject *dreampy_get_other_info(
    PyObject* /*self*/, PyObject *args, PyObject *kwargs
) {
    static const char *kwlist[] = {"ptr", "name", NULL};
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

    DREAM::OtherQuantity *o = get_other(sim, name);

    if (o == NULL)
        return NULL;

    PyObject *dct = PyDict_New();

    PyDict_SetItemString(dct, "description", PyUnicode_FromString(o->GetDescription().c_str()));
    PyDict_SetItemString(dct, "nelements", PyLong_FromUnsignedLongLong(o->NumberOfElements()));
    PyDict_SetItemString(dct, "nmultiples", PyLong_FromUnsignedLongLong(o->NumberOfMultiples()));

    return dct;
}

}
