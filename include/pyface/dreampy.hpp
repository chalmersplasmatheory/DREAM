#ifndef _DREAM_PYFACE_HPP
#define _DREAM_PYFACE_HPP

#define PY_SSIZE_T_CLEAN
#include <Python.h>

extern "C" {
    static PyObject *dreampy_run(PyObject*, PyObject*);
}

#endif/*_DREAM_PYFACE_HPP*/
