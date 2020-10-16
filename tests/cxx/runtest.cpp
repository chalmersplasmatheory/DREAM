/* Unit tests */

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <petsc.h>
#include <fenv.h>

#include <softlib/SOFTLibException.h>

#include "tests/cxx/config.h"
#include "FVM/FVMException.hpp"

#include "UnitTest.hpp"

// Tests
#include "tests/DREAM/BoundaryFlux.hpp"
#include "tests/DREAM/IonRateEquation.hpp"
#include "tests/DREAM/RunawayFluid.hpp"
#include "tests/DREAM/AvalancheSourceRP.hpp"
#include "tests/DREAM/MeanExcitationEnergy.hpp"

#include "tests/FVM/AdvectionTerm.hpp"
#include "tests/FVM/AdvectionDiffusionTerm.hpp"
#include "tests/FVM/AnalyticBRadialGridGenerator.hpp"
#include "tests/FVM/DiffusionTerm.hpp"
#include "tests/FVM/Grid.hpp"
#include "tests/FVM/Interpolator1D.hpp"
#include "tests/FVM/Interpolator3D.hpp"
#include "tests/FVM/PXiExternalKineticKinetic.hpp"

using namespace std;
using namespace DREAMTESTS;

vector<UnitTest*> tests;

void add_test(UnitTest *t) {
	tests.push_back(t);
}
void init() {
    add_test(new DREAMTESTS::_DREAM::AvalancheSourceRP("dream/avalanche"));
    add_test(new DREAMTESTS::_DREAM::BoundaryFlux("dream/boundaryflux"));
    add_test(new DREAMTESTS::_DREAM::IonRateEquation("dream/ionrateequation"));
    add_test(new DREAMTESTS::_DREAM::MeanExcitationEnergy("dream/meanexcitationenergy"));
    add_test(new DREAMTESTS::_DREAM::RunawayFluid("dream/runawayfluid"));

    add_test(new DREAMTESTS::FVM::AdvectionTerm("fvm/advectionterm"));
    add_test(new DREAMTESTS::FVM::AdvectionDiffusionTerm("fvm/advectiondiffusionterm"));
    add_test(new DREAMTESTS::FVM::AnalyticBRadialGridGenerator("fvm/fluxsurfaceaverage"));
    add_test(new DREAMTESTS::FVM::DiffusionTerm("fvm/diffusionterm"));
    add_test(new DREAMTESTS::FVM::Grid("fvm/grid"));
    add_test(new DREAMTESTS::FVM::Interpolator1D("fvm/interpolator1d"));
    add_test(new DREAMTESTS::FVM::Interpolator3D("fvm/interpolator3d"));
    add_test(new DREAMTESTS::FVM::PXiExternalKineticKinetic("fvm/boundaryflux/2kinetic"));
}

/**
 * Check if a test with the given name
 * exists in the test list "tests".
 *
 * name: Name of test to look for.
 *
 * RETURNS the index of test if found,
 * otherwise, -1.
 */
int has_test(const char *name) {
	size_t i;
	for (i = 0; i < tests.size(); i++) {
		if (tests[i]->HasName(name))
			return i;
	}

	return -1;
}

/**
 * Run the test with index 'index' in
 * the test list "tests".
 *
 * index: Index in "tests" of the test
 *    to run.
 * runningAll: TRUE if all (or a group of)
 *    tests are being run.
 *
 * RETURNS 1 if the test failed, 0 otherwise.
 */
int run_test(int index, bool runningAll) {
    cout << "\x1B[1m:: " << tests[index]->GetName() << "\x1B[0m" << endl;
    try {
        if (tests[index]->Run(runningAll)) {
            cout << "\x1B[1;32m[SUCCESS]\x1B[0m Test '" << tests[index]->GetName() << "' completed successfully." << endl;
            return 0;
        } else {
            cout << "\x1B[1;31m[FAIL]\x1B[0m    Test '" << tests[index]->GetName() << "' failed." << endl;
            return 1;
        }
    } catch (DREAM::FVM::FVMException& ex) {
        cout << "\x1B[1;31m[ERROR]\x1B[0m   --> " << ex.whats() << endl;
        cout << "\x1B[1;31m[FAIL]\x1B[0m    Test '" << tests[index]->GetName() << "' unexpectedly threw exception." << endl;
    } catch (SOFTLibException& ex) {
        cout << "\x1B[1;31m[ERROR]\x1B[0m   --> " << ex.whats() << endl;
        cout << "\x1B[1;31m[FAIL]\x1B[0m    Test '" << tests[index]->GetName() << "' unexpectedly threw exception." << endl;
    }

    return 1;
}
/**
 * Runs all tests in the list "tests".
 *
 * RETURNS the number of tests that failed.
 */
int run_all_tests() {
	size_t i;
	int failed = 0;

	for (i = 0; i < tests.size(); i++) {
		failed += run_test(i, true);
	}

	return failed;
}

/**
 * Print help text.
 */
void help() {
    printf(
        "Unit tests for DREAM\n\n"

        "Usage:\n"
        "    dreamtests          Show this help message.\n"
        "    dreamtests all      Run all tests.\n"
        "    dreamtests [test1 [test2 [...]]]\n"
        "                        Runs the tests with names 'test1', 'test2' etc.\n\n"

        "Available tests:\n"
    );
    for (unsigned int i = 0; i < tests.size(); i++) {
        printf("    %s\n", tests[i]->GetName().c_str());
    }

    printf("\n");
}

int main(int argc, char *argv[]) {
	int i, t, failed=0, total=0;

    PetscInitialize(&argc, &argv, NULL, NULL);

	init();

    // Enable floating-point exceptions
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

	if (argc == 1) {
        help();
    } else if (argc == 2 && !strcmp(argv[1], "all")) {
		failed = run_all_tests();
        total = tests.size();
	} else {
        total = argc-1;
		for (i = 1, failed=0; i < argc; i++) {
			if ((t=has_test(argv[i]))>=0) {
				failed += run_test(t, (argc>2));
			} else {
				cerr << "\x1B[1;31m[ERROR]\x1B[0m   Unrecognized test: '" << argv[i] << "'." << endl;
				failed++;
			}
		}
	}

	cout << (total-failed) << " of " << total << " tests passed." << endl;

    PetscFinalize();

	return failed;
}

