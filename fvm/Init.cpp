/**
 * Initialization for the FVM library.
 */

#include <gsl/gsl_errno.h>
#include "FVM/GSLException.hpp"
#include "FVM/Init.hpp"

void dream_fvm_initialize() {
    dream_fvm_set_gsl_error_handler();
}
void dream_fvm_finalize() {}

/**
 * Set up a GSL error handler which converts a GSL
 * error into a C++ exception.
 */
void dream_fvm_set_gsl_error_handler() {
    gsl_set_error_handler(&dream_fvm_gsl_error_handler);
}

/**
 * Definition of the DREAM GSL error handler which
 * converts a GSL error into a C++ exception.
 *
 * reason:    Descriptive error message from GSL.
 * file:      Name of file in which the error occured.
 * line:      Name of line in 'file' on which the error occured.
 * gsl_errno: Internal code used by GSL to identify the error.
 */
void dream_fvm_gsl_error_handler(
    const char *reason, const char *file, int line, int gsl_errno
) {
    throw DREAM::FVM::GSLException(reason, file, line, gsl_errno);
}

