#ifndef _DREAM_NUMPY_H
#define _DREAM_NUMPY_H

/**
 * Import NumPy properly.
 */

#define PY_ARRAY_UNIQUE_SYMBOL DREAM_PyArray_API

// Macro for translation unit
#ifndef INIT_NUMPY_ARRAY_CPP
#   define NO_IMPORT_ARRAY
#endif

#ifndef NPY_NO_DEPRECATED_API
#   define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

// Include numpy
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

#endif/*_DREAM_NUMPY_H*/
