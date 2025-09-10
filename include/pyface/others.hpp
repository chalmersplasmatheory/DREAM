#ifndef _DREAM_PYFACE_OTHERS_HPP
#define _DREAM_PYFACE_OTHERS_HPP

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>
#include "DREAM/OtherQuantity.hpp"
#include "DREAM/Simulation.hpp"

DREAM::OtherQuantity *get_other(DREAM::Simulation*, const char*);

extern "C" {
    static PyObject *dreampy_get_others(PyObject*, PyObject*);
    static PyObject *dreampy_get_other_data(PyObject*, PyObject*, PyObject*);
    static PyObject *dreampy_get_other_info(PyObject*, PyObject*, PyObject*);
}

#endif/*_DREAM_PYFACE_OTHERS_HPP*/
