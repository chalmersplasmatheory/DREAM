#ifndef _DREAM_PYFACE_HPP
#define _DREAM_PYFACE_HPP

#define PY_SSIZE_T_CLEAN
#include <Python.h>

extern "C" {
    static PyObject *dreampy_run(PyObject*, PyObject*);
    static PyObject *dreampy_run_simulation(PyObject*, PyObject*);
    static PyObject *dreampy_setup_simulation(PyObject*, PyObject*);
    static PyObject *dreampy_register_callback_timestep_finished(PyObject*, PyObject*, PyObject*);
    static PyObject *dreampy_register_callback_iteration_finished(PyObject*, PyObject*, PyObject*);
}

#endif/*_DREAM_PYFACE_HPP*/
