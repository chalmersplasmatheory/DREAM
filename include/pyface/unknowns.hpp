#ifndef _DREAM_PYFACE_UNKNOWNS_HPP
#define _DREAM_PYFACE_UNKNOWNS_HPP

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>

extern "C" {
    static PyObject *dreampy_get_unknowns(PyObject*, PyObject*);
    static PyObject *dreampy_get_unknown_info(PyObject*, PyObject*, PyObject*);
}

#endif/*_DREAM_PYFACE_UNKNOWNS_HPP*/
