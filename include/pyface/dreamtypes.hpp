#ifndef _DREAM_PYFACE_TYPES_HPP
#define _DREAM_PYFACE_TYPES_HPP

#include "DREAM/config.h"

/**
 * Template function for converting array data
 * from a type 'T1' to a type 'T2'.
 *
 * inp:  Input data.
 * ndim: Number of dimensions of data.
 * dims: Array holding size of dimensions.
 */
template<typename T1, typename T2>
T1 *dreampy_convert(T2 *inp, int ndim, npy_intp *dims, npy_intp *strides) {
    // Calculate total size of array
    len_t size = 1;
    for (int i = 0; i < ndim; i++)
        size *= dims[i];
	len_t stride = strides[ndim-1] / sizeof(T2);

    T1 *out = new T1[size];
    for (len_t i = 0; i < size; i++)
        out[i] = inp[i*stride];

    return out;
}

#endif/*_DREAMPY_PYFACE_TYPES_HPP*/

