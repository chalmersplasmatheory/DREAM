#ifndef _DREAM_PYFACE_SIMULATION_HPP
#define _DREAM_PYFACE_SIMULATION_HPP

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>
#include "DREAM/Simulation.hpp"

DREAM::Simulation *get_simulation_from_capsule(PyObject*);

extern "C" {
    static PyObject *dreampy_get_current_time(PyObject*, PyObject*);
    static PyObject *dreampy_get_max_time(PyObject*, PyObject*);
    static PyObject *dreampy_get_radius_vector(PyObject*, PyObject*);
    static PyObject *dreampy_get_time_vector(PyObject*, PyObject*);
}

#endif/*_DREAM_PYFACE_SIMULATION_HPP*/
