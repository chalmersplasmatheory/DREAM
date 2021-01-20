/**
 * Implementation of tests for the combined advection & diffusion term.
 */

#include "FVM/Equation/AdvectionDiffusionTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Equation/BoundaryConditions/PXiInternalTrapping.hpp"
#include "FVM/Equation/Operator.hpp"
#include "AdvectionDiffusionTerm.hpp"
#include "GeneralAdvectionTerm.hpp"
#include "GeneralDiffusionTerm.hpp"
#include "GeneralAdvectionDiffusionTerm.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Check if the implementation of the combined advection
 * & diffusion term preserves density.
 */
bool AdvectionDiffusionTerm::CheckConservativity(DREAM::FVM::Grid *grid) {
    bool isConservative = true;

    DREAM::FVM::Operator *Op = new DREAM::FVM::Operator(grid);
    Op->AddTerm(new GeneralAdvectionTerm(grid,1.0));
    Op->AddTerm(new GeneralDiffusionTerm(grid,1.0));
    if(grid->HasTrapped())
        Op->AddBoundaryCondition(new DREAM::FVM::BC::PXiInternalTrapping(grid,Op));
    Op->SetAdvectionBoundaryConditions(
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET, 
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET
    );

    const len_t ncells = grid->GetNCells();
    const len_t NNZ_PER_ROW = Op->GetNumberOfNonZerosPerRow();
    DREAM::FVM::Matrix *mat = new DREAM::FVM::Matrix(ncells, ncells, NNZ_PER_ROW);

    /**
     * Successively build
     *
     *   i = 0: Only r-terms
     *   i = 1: Only p1-terms
     *   i = 2: Only p2-terms
     *   i = 3: All terms
     */
    for (len_t i = 0; i < 6; i++) {
        Op->RebuildTerms(5-i, 0, nullptr);
        Op->SetMatrixElements(mat, nullptr);
        mat->Assemble();

        const real_t TOLERANCE = 2*NNZ_PER_ROW*ncells * std::numeric_limits<real_t>::epsilon();

        if (!IsConservative(mat, grid, TOLERANCE)) {
            const char *dim = (i==0?"r" : (i==1?"p1" : (i==2?"p2":"every")));
            this->PrintError("Combined advection-diffusion term is not conservative in '%s' component.", dim);

            isConservative = false;
        }
        mat->Zero();
    }

    delete mat;
    delete Op;

    return isConservative;
}

/**
 * Check that this term is evaluated correctly
 * (we only override it, but don't actually implement it)
 */
bool AdvectionDiffusionTerm::CheckValue(DREAM::FVM::Grid*) {
    return false;
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

