/**
 * Implementation of tests for the 'DiffusionTerm' class
 * in the DREAM FVM library.
 */

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Matrix.hpp"
#include "DiffusionTerm.hpp"
#include "GeneralDiffusionTerm.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Check if the implementation of the diffusion term
 * preserves density.
 */
bool DiffusionTerm::CheckConservativity(DREAM::FVM::Grid *grid) {
    bool isConservative = true;
    GeneralDiffusionTerm *gdt = new GeneralDiffusionTerm(grid);

    const len_t ncells = grid->GetNCells();
    const len_t NNZ_PER_ROW = 11;
    DREAM::FVM::Matrix *mat = new DREAM::FVM::Matrix(ncells, ncells, NNZ_PER_ROW);

    for (len_t i = 0; i < 5; i++) {
        gdt->Rebuild(i);
        gdt->SetMatrixElements(mat, nullptr);
        mat->Assemble();

        const real_t TOLERANCE = NNZ_PER_ROW*ncells * std::numeric_limits<real_t>::epsilon();

        if (!IsConservative(mat, grid, TOLERANCE)) {
            const char *dim = (i==0?"rr" : (i==1?"11" : (i==2?"22": (i==3?"12" : (i==4?"21":"every")))));
            this->PrintError("Diffusion term is not conservative in '%s' component.", dim);

            isConservative = false;
        }

        mat->Zero();
    }

    delete mat;
    delete gdt;

    return isConservative;
}

/**
 * Run all tests for this module.
 */
bool DiffusionTerm::Run(bool) {
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
        this->PrintOK("The general diffusion term conserves density");

    return success;
}

