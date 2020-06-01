
#define INIT_NUMPY_ARRAY_CPP
#include "pyface/numpy.h"

/**
 * Initialize numpy (must (allegedly) be done in
 * every _file_ using the NumPy C API)
 */
int init_numpy() {
    import_array();
    return 0;
}

const static int numpy_initialized = init_numpy();

