/**
 * Initialization routines for the TQS library.
 */

#include <petsc.h>
#include "TQS/Init.h"


/**
 * Initializes the TQS library.
 */
void tqs_initialize(int *argc, char **argv[]) {
    PetscInitialize(argc, argv, NULL, NULL);
}

/**
 * De-initialize TQS and release all used resources.
 */
void tqs_finalize() {
    PetscFinalize();
}

