/**
 * Python interface to DREAM.
 */

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>
#include <iostream>
#include "DREAM/Init.h"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "pyface/dreampy.hpp"
#include "pyface/settings.hpp"


static PyMethodDef dreampyMethods[] = {
    {"run", dreampy_run, METH_VARARGS, "Run DREAM."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef dreampyModule = {
    PyModuleDef_HEAD_INIT,
    "libdreampy",           // Name of module
    nullptr,                // Module documentation (may be NULL)
    -1,                     // Size of per-interpreter state of the module,
                            // of -1 if the module keeps state in global variables
    dreampyMethods,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/**
 * Module initialization.
 */
PyMODINIT_FUNC PyInit_libdreampy() {
    return PyModule_Create(&dreampyModule);
}

/**
 * Run a DREAM simulation. This function takes a Python
 * dictionary with the simulation settings as input.
 */
extern "C" {
static PyObject *dreampy_run(PyObject* /*self*/, PyObject *args) {
    bool success = true;

    PyObject *dict;
    if (!PyArg_ParseTuple(args, "O", &dict))
        return NULL;

    // Initialize the DREAM kernel
    dream_initialize();

    try {
        DREAM::Settings *settings = dreampy_loadsettings(dict);
        //DREAM::Simulation *sim = DREAM::SimulationGenerator::ProcessSettings(settings);
        //sim->Run();

        // TODO save output
    } catch (DREAM::FVM::FVMException& ex) {
        PyErr_SetString(PyExc_RuntimeError, ex.what());
        success = false;
    }

    // De-initialize the DREAM kernel
    dream_finalize();

    /*printf("Step 1\n");
    PyObject *keys = PyDict_Keys(dict);
    Py_ssize_t len = PyList_Size(keys);

    printf("Keys:\n");
    for (Py_ssize_t i = 0; i < len; i++) {
        PyObject *obj = PyList_GetItem(keys, i);
        Py_ssize_t slen;
        const char *s = PyUnicode_AsUTF8AndSize(obj, &slen);

        std::cout << "  " << s << std::endl;
    }*/

    if (success)
        Py_RETURN_NONE;
    else
        return NULL;
}
}

