/**
 * Implementation of UnitTest class.
 */

#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <string>
#include <softlib/SFile.h>

#include "DREAM/config.h"
#include "UnitTest.hpp"

#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"

using namespace std;
using namespace DREAMTESTS;

/**
 * Constructor
 */
UnitTest::UnitTest(const string& name) {
    this->name = name;
    SeedRand();
}

/**
 * Get name of this unit test.
 */
string& UnitTest::GetName() { return this->name; }

/**
 * Check if this unit test has
 * the given name.
 *
 * cmp: Name to check against.
 */
bool UnitTest::HasName(const string& cmp) { return (this->name==cmp); }

/**
 * Print an [ERROR] message.
 * Uses printf syntax.
 */
void UnitTest::PrintError(const string s, ...) {
    va_list args;
    va_start(args, s);

    fprintf(stderr, "\x1B[1;31m[ERROR]\x1B[0m   ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
/**
 * Print an [OK] message.
 * Uses printf syntax.
 */
void UnitTest::PrintOK(const string s, ...) {
    va_list args;
    va_start(args, s);

    fprintf(stderr, "\x1B[1;32m[OK]\x1B[0m      --> ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
/**
 * Print a general status message.
 * Uses printf syntax.
 */
void UnitTest::PrintStatus(const string s, ...) {
    va_list args;
    va_start(args, s);

	fprintf(stderr, "\x1B[1m[INFO]\x1B[0m    --> ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
/**
 * Print a [WARNING] message.
 * Uses printf syntax.
 */
void UnitTest::PrintWarning(const string s, ...) {
    va_list args;
    va_start(args, s);

    fprintf(stderr, "\x1B[1;93m[WARNING]\x1B[0m ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}

/**
 * Generate a random number between
 * 0.0 and 1.0 (inclusive in both limits).
 */
real_t UnitTest::Rand() {
	return (((real_t)rand()) / ((real_t)RAND_MAX));
}

/**
 * Seed the random number generator.
 */
void UnitTest::SeedRand() {
    srand(time(NULL));
}

/**
 * Save the given 1D array to the named file. The variable
 * is named "data" in the output file.
 *
 * filename: Name of file to save variable to.
 * var:      1D array to save.
 * n:        Number of elements in 'var'.
 */
void UnitTest::SaveVariable(const string &filename, real_t *var, const len_t n) {
	SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
	sf->WriteList("data", var, n);
	sf->Close();

	delete sf;
}

/**
 * Save the given 2D array to the named file. The variable
 * is named "data" in the output file.
 *
 * filename: Name of file to save variable to.
 * var:      2D array to save.
 * m:        Number of rows in 'var'.
 * n:        Number of columns in 'var'.
 */
void UnitTest::SaveVariable(const string &filename, real_t **var, const len_t m, const len_t n) {
	SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
	sf->WriteArray("data", var, m, n);
	sf->Close();

	delete sf;
}

