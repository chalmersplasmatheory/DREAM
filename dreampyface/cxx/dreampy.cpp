/**
 * Python interface to DREAM.
 */
#include <csignal>

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>
#include <iostream>
#include <vector>
#include "DREAM/Init.h"
#include "DREAM/OutputGeneratorSFile.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "pyface/callback.hpp"
#include "pyface/dreampy.hpp"
#include "pyface/others.hpp"
#include "pyface/settings.hpp"
#include "pyface/simulation.hpp"
#include "pyface/SFile_Python.hpp"
#include "pyface/unknowns.hpp"
#include "softlib/Timer.h"
#include "DREAM/IO.hpp"
#include "DREAM/QuitException.hpp"


static PyMethodDef dreampyMethods[] = {
    {"get_current_time", dreampy_get_current_time, METH_VARARGS, "Returns the current time of the given simulation."},
    {"get_max_time", dreampy_get_max_time, METH_VARARGS, "Returns the maximum simulation time of the given simulation."},
    {"get_others", dreampy_get_others, METH_VARARGS, "Returns a dictionary with information about all other quantities of the equation system."},
    {"get_other_data", (PyCFunction)(void(*)(void))dreampy_get_other_data, METH_VARARGS | METH_KEYWORDS, "Returns a dictionary with data for the other quantity corresponding to all completed time steps."},
    {"get_other_info", (PyCFunction)(void(*)(void))dreampy_get_other_info, METH_VARARGS | METH_KEYWORDS, "Returns a dictionary with information about the named other quantity."},
    {"get_radius_vector", dreampy_get_radius_vector, METH_VARARGS, "Returns the radial grid used in the simulation."},
    {"get_time_vector", dreampy_get_time_vector, METH_VARARGS, "Returns a list of time steps completed so far."},
    {"get_unknowns", dreampy_get_unknowns, METH_VARARGS, "Returns a dictionary with information about all unknowns of the equation system."},
    {"get_unknown_data", (PyCFunction)(void(*)(void))dreampy_get_unknown_data, METH_VARARGS | METH_KEYWORDS, "Returns a dictionary with data for the unknown quantity corresponding to all completed time steps."},
    {"get_unknown_info", (PyCFunction)(void(*)(void))dreampy_get_unknown_info, METH_VARARGS | METH_KEYWORDS, "Returns a dictionary with information about the named unknown quantity."},
    {"register_callback_timestep_finished", (PyCFunction)(void(*)(void))dreampy_register_callback_timestep_finished, METH_VARARGS | METH_KEYWORDS, "Register a function to call when a timestep has been completed."},
    {"register_callback_iteration_finished", (PyCFunction)(void(*)(void))dreampy_register_callback_iteration_finished, METH_VARARGS | METH_KEYWORDS, "Register a function to call when a solver iteration has finished."},
    {"run", dreampy_run, METH_VARARGS, "Run DREAM."},
    {"run_simulation", dreampy_run_simulation, METH_VARARGS, "Run a previously constructed DREAM simulation object."},
    {"setup_simulation", dreampy_setup_simulation, METH_VARARGS, "Construct a DREAM simulation object."},
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
static PyObject *dreampy_run(PyObject *self, PyObject *args) {
    PyObject *sim = dreampy_setup_simulation(self, args);
    if (sim == NULL)
        return NULL;

    PyObject *_args = PyTuple_Pack(1, sim);
    PyObject *out = dreampy_run_simulation(self, _args);
    Py_DECREF(_args);

    return out;
    /*bool success = true;

    PyObject *dict;
    if (!PyArg_ParseTuple(args, "O", &dict)) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "The argument to 'run' must be a Python dictionary."
        );
        return NULL;
    }

    // Initialize the DREAM kernel
    dream_initialize();

    DREAM::OutputGeneratorSFile *ogs;
    SFile_Python *sfp = new SFile_Python();
    DREAM::Settings *settings;
    DREAM::Simulation *sim;

    try {
        settings = dreampy_loadsettings(dict);
        sim = DREAM::SimulationGenerator::ProcessSettings(settings);
        
        ogs = new DREAM::OutputGeneratorSFile(
            sim->GetEquationSystem(), sfp
        );
        sim->SetOutputGenerator(ogs);

        // Register callback functions
        register_callback_functions(sim);

        sim->Run();

        // Generate output
        Timer t;
        sim->Save();
        std::cout << "Time to save output: " << t.ToString() << std::endl;
    } catch (DREAM::FVM::FVMException& ex) {
        PyErr_SetString(PyExc_RuntimeError, ex.what());
        success = false;
    }

    // De-initialize the DREAM kernel
    dream_finalize();

    if (success) {
        PyObject *d = sfp->GetPythonDict();
        
        // TODO delete simulation

        return d;
    } else
        return NULL;*/
}

/**
 * Set up a simulation object from the given settings
 * (passed as a Python dictionary).
 */
static PyObject *dreampy_setup_simulation(
    PyObject* /*self*/, PyObject *args
) {
    bool success = true;

    PyObject *dict;
    if (!PyArg_ParseTuple(args, "O", &dict)) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "The argument to 'run' must be a Python dictionary."
        );
        return NULL;
    }

    // Initialize the DREAM kernel
    dream_initialize();

    DREAM::Settings *settings;
    DREAM::Simulation *sim;

    try {
        settings = dreampy_loadsettings(dict);
        sim = DREAM::SimulationGenerator::ProcessSettings(settings);
        
        
        sim->GetEquationSystem()->GetTimeStepper()->SetPythonCaller(dreampy_callback_return_bool);
    } catch (DREAM::FVM::FVMException& ex) {
        PyErr_SetString(PyExc_RuntimeError, ex.what());
        success = false;
    }
    
    if (success)
        return PyCapsule_New(sim, "sim", NULL);
    else
        return NULL;
}

/**
 * Handle a 'SIGQUIT' signal.
 */
void sig_quit(int) {
    throw DREAM::QuitException("The user requested execution to stop.");
}

/**
 * Handle a 'SIGQUIT' signal.
 */
void sig_segv(int) {
    throw DREAM::QuitException("A segmentation fault occured.");
}

/**
 * Run a previously constructed Simulation, passed
 * as a PyCapsule object to this function.
 */
static PyObject *dreampy_run_simulation(
    PyObject* /*self*/, PyObject *args
) {
    dream_make_splash();
    bool success = true;
    DREAM::Simulation *sim = get_simulation_from_capsule(args);
    DREAM::OutputGeneratorSFile *ogs = nullptr;
    SFile_Python *sfp = new SFile_Python();

    // Register callback functions
    register_callback_functions(sim);

    // Allow the user to press Ctrl+\ or Ctrl+Y to quit the simulation early
    PetscPopSignalHandler();
    std::signal(SIGQUIT, sig_quit);
    std::signal(SIGSEGV, sig_segv);

    try {
        ogs = new DREAM::OutputGeneratorSFile(
            sim->GetEquationSystem(), sfp
        );
        sim->Run();

        // Generate output
        Timer t;
        ogs->Save();
        std::cout << "Time to save output: " << t.ToString() << std::endl;
    } catch (DREAM::FVM::FVMException& ex) {
        PyErr_SetString(PyExc_RuntimeError, ex.what());
        success = false;
    }
    
    if (sim != nullptr) {
        try {
            sim->Save();
        } catch (H5::FileIException &ex) {
            PyErr_SetString(PyExc_RuntimeError, ex.getDetailMsg().c_str());
        }
    }

    // De-initialize the DREAM kernel
    dream_finalize();

    if (success) {
        SFile_Python *sfp = static_cast<SFile_Python*>(ogs->GetSFile());
        PyObject *d = sfp->GetPythonDict();
        
        // TODO delete simulation

		if (ogs != nullptr)
			delete ogs;

        delete sim;
        return d;
    } else {
		if (ogs != nullptr)
			delete ogs;

        delete sim;
        return NULL;
    }
}

/**
 * Python signature:
 *
 *   def register_callback_timestep_finished(func)
 *
 * where
 *
 *   func: Callback function to register.
 */
static PyObject *dreampy_register_callback_timestep_finished(
    PyObject* /*self*/, PyObject *args, PyObject *kwargs
) {
    static const char *kwlist[] = {"func", NULL};
    PyObject *func;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", const_cast<char**>(kwlist), &func))
        return NULL;

    // Verify that 'func' is indeed a function
    if (!PyFunction_Check(func)) {
        // Raise exception
        PyErr_SetString(
            PyExc_RuntimeError,
            "The argument 'func' must be a Python function."
        );
        return NULL;
    }

    // Register function
    Py_INCREF(func);
    callback_timestep.push_back(func);

    Py_RETURN_NONE;
}

static PyObject *dreampy_register_callback_iteration_finished(
    PyObject* /*self*/, PyObject *args, PyObject *kwargs
) {
    static const char *kwlist[] = {"func", NULL};
    PyObject *func;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", const_cast<char**>(kwlist), &func))
        return NULL;

    // Verify that 'func' is indeed a function
    if (!PyFunction_Check(func)) {
        // Raise exception
        PyErr_SetString(
            PyExc_RuntimeError,
            "The argument 'func' must be a Python function."
        );
        return NULL;
    }

    // Register function
    Py_INCREF(func);
    callback_iteration.push_back(func);

    Py_RETURN_NONE;
}

}//extern "C"

#include "others.cpp"
#include "simulation.cpp"
#include "unknowns.cpp"

