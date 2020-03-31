/**
 * Implementation of tests for the 'AdvectionTerm' class
 * in the DREAM FVM library.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Matrix.hpp"
#include "AdvectionTerm.hpp"
#include "GeneralAdvectionTerm.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Check the implementation of the advection term
 * preserves density.
 */
bool AdvectionTerm::CheckConservativity(DREAM::FVM::Grid *grid) {
    bool isConservative = true;
    GeneralAdvectionTerm *gat = new GeneralAdvectionTerm(grid);

    const len_t ncells = grid->GetNCells();
    const len_t NNZ_PER_ROW = 7;
    DREAM::FVM::Matrix *mat = new DREAM::FVM::Matrix(ncells, ncells, NNZ_PER_ROW);

    for (len_t i = 0; i < 3; i++) {
        gat->Rebuild(i);
        gat->SetMatrixElements(mat);

        mat->Assemble();

        const real_t TOLERANCE = NNZ_PER_ROW*ncells * std::numeric_limits<real_t>::epsilon();

        if (!IsConservative(mat, grid, TOLERANCE)) {
            const char *dim = (i==0?"r" : (i==1?"p1" : (i==2?"p2":"every")));
            this->PrintError("Advection term is not conservative in '%s' dimension.", dim);

            isConservative = false;
        }

        mat->Zero();
    }

    delete mat;
    delete gat;

    return isConservative;
}

/**
 * Run all tests for this module.
 */
bool AdvectionTerm::Run(bool) {
    bool success = true;

    // Check the conservativity of the operator on
    // ALL implemented grid combinations
    // (see 'EquationTerm.cpp' for implementation of this
    // routine; see 'UnitTest.cpp' for implementation of
    // the various grid combinations)
    if (!this->EquationTerm::CheckConservativity()) {
        this->PrintError("Conservativity test failed");
        success = false;
    } else
        this->PrintOK("The general advection term conserves density");

    return success;
}

