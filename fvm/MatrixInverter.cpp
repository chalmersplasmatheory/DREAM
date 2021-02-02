/**
 * Common routines for the 'MatrixInverter' classes.
 */

#include <iostream>
#include "FVM/MatrixInverter.hpp"


using namespace DREAM::FVM;
using namespace std;


/**
 * Print info about the most recently factored matrix.
 */
void MatrixInverter::PrintInfo() {
    Mat mat;
    MatInfo info;
    PC pc;

    KSPGetPC(this->ksp, &pc);
    PCFactorGetMatrix(pc, &mat);
    MatGetInfo(mat, MAT_GLOBAL_MAX, &info);

    cout << ":: FACTORED MATRIX INFORMATION" << endl;
    cout << "   Block size:             " << info.block_size << endl;
    cout << "   Non-zeros" << endl;
    cout << "     allocated:            " << info.nz_allocated << endl;
    cout << "     used:                 " << info.nz_used << endl;
    cout << "     unneeded:             " << info.nz_unneeded << endl;
    cout << "   Memory allocated:       " << info.memory << endl;
    cout << "   Matrix assemblies:      " << info.assemblies << endl;
    cout << "   Number of mallocs:      " << info.mallocs << endl;
    cout << "   Fill ratio for LU/ILU" << endl;
    cout << "     given:                " << info.fill_ratio_given << endl;
    cout << "     needed:               " << info.fill_ratio_needed << endl;
    cout << "   # mallocs during fact.: " << info.factor_mallocs << endl;
}

