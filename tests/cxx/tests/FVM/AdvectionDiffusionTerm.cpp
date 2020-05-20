/**
 * Implementation of tests for the combined advection & diffusion term.
 */

#include "FVM/Equation/AdvectionDiffusionTerm.hpp"
#include "FVM/Matrix.hpp"
#include "AdvectionDiffusionTerm.hpp"
#include "GeneralAdvectionDiffusionTerm.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Check if the implementation of the combined advection
 * & diffusion term preserves density.
 */
bool AdvectionDiffusionTerm::CheckConservativity(DREAM::FVM::Grid *grid) {
    bool isConservative = true;
    GeneralAdvectionDiffusionTerm *gadt = new GeneralAdvectionDiffusionTerm(grid);

    const len_t ncells = grid->GetNCells();
    const len_t NNZ_PER_ROW = gadt->GetNumberOfNonZerosPerRow();
    DREAM::FVM::Matrix *mat = new DREAM::FVM::Matrix(ncells, ncells, NNZ_PER_ROW);

    /**
     * Successively build
     *
     *   i = 0: Only r-terms
     *   i = 1: Only p1-terms
     *   i = 2: Only p2-terms
     *   i = 3: All terms
     */
    for (len_t i = 0; i < 4; i++) {
        gadt->Rebuild(i, 0, nullptr);
        gadt->SetMatrixElements(mat, nullptr);
        mat->Assemble();

        const real_t TOLERANCE = 2*NNZ_PER_ROW*ncells * std::numeric_limits<real_t>::epsilon();

        if (!IsConservative(mat, grid, TOLERANCE)) {
            const char *dim = (i==0?"r" : (i==1?"p1" : (i==2?"p2":"every")));
            this->PrintError("Combined advection-diffusion term is not conservative in '%s' component.", dim);

            isConservative = false;
        }

        if (i == 1)
            mat->View(DREAM::FVM::Matrix::BINARY_MATLAB);

        mat->Zero();
    }

    delete mat;
    delete gadt;

    return isConservative;
}

/**
 * Run all tests for this module.
 */
bool AdvectionDiffusionTerm::Run(bool) {
    bool success = true;

    if (!this->EquationTerm::CheckConservativity()) {
        this->PrintError("Conservativity test failed");
        success = false;
    } else
        this->PrintOK("The general combined advection-diffusion term conserves density.");

    return success;
}

