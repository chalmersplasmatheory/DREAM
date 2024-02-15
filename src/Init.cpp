/**
 * Initialization routines for the DREAM library.
 */

#include <petsc.h>
#include "DREAM/Init.h"
#include "FVM/Init.hpp"
#include <iostream>
using namespace std;
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

void dream_make_splash(){
    cout << endl;
    cout << R"( It's time to...)" << endl;
    cout << endl;
    cout << R"( * * ________   _____  ______ ___* *  ___  ___       )" << endl;
    cout << R"(  * \\   __   \/ __  \/  ____//   \ //   \/   \      )" << endl;
    cout << R"(    //  / //  / /_// /  /__ // /\  \/  / / /  /   *  )" << endl;
    cout << R"(   //  / //  /      /   __/// /_/  /  / / /  /   * * )" << endl;
    cout << R"( _//  /_//  /  /\  \   /__//  __  /  / / /  /     *  )" << endl;
    cout << R"( \\________/__/ \\__\_____//_///_/__//_//__/     *   )" << endl;
    cout << R"(               * *          * *      *         *     )" << endl;
    cout << R"(      * *       *       * *      * *             *   )" << endl;
    cout << R"(     * * *                                 ...baby...)" << endl;
    cout << endl;
}