/**
 * Initialization routines for the DREAM library.
 */

#include <petsc.h>
#include "DREAM/Init.h"
#include "FVM/Init.hpp"


/**
 * Initializes the DREAM library.
 */
void dream_initialize() {
    dream_initialize(nullptr, nullptr);
}
void dream_initialize(int *argc, char **argv[]) {
    if (argc == nullptr)
        PetscInitializeNoArguments();
    else
        PetscInitialize(argc, argv, NULL, NULL);

    dream_fvm_initialize();

    //PetscInfoAllow(PETSC_TRUE);
    //PetscInfoSetFile("petsc_out.txt", "w");
}

/**
 * De-initialize DREAM and release all used resources.
 */
void dream_finalize() {
    dream_fvm_finalize();
    PetscFinalize();
}

