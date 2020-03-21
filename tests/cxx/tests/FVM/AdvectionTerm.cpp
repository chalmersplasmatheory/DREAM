/**
 * Implementation of tests for the 'AdvectionTerm' class
 * in the TQS FVM library.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Matrix.hpp"
#include "AdvectionTerm.hpp"
#include "GeneralAdvectionTerm.hpp"


using namespace TQSTESTS::FVM;

/**
 * Check the implementation of the advection term
 * preserves density.
 */
bool AdvectionTerm::CheckConservativity(TQS::FVM::RadialGrid *rg) {
    GeneralAdvectionTerm *gat = new GeneralAdvectionTerm(rg);

    const len_t ncells = rg->GetNCells();
    TQS::FVM::Matrix *mat = new TQS::FVM::Matrix(ncells, ncells, 3);

    gat->Rebuild(0);
    gat->SetMatrixElements(mat);

    mat->Assemble();

    const real_t TOLERANCE = 3*ncells * std::numeric_limits<real_t>::epsilon();
    bool isConservative = IsConservative(mat, rg, TOLERANCE);

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

