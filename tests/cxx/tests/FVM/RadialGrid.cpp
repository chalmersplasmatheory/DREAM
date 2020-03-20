/**
 * Test for the 'RadialGrid' in the FVM library.
 */

#include "UnitTest.hpp"
#include "RadialGrid.hpp"


using namespace TQSTESTS::FVM;

/**
 * Check whether a general radial grid can be appropriately
 * constructed.
 */
bool RadialGrid::CheckGeneralGrid() {
    bool success = true;

    TQS::FVM::RadialGrid *rg = this->InitializeGeneralGridPXi();

    ///////////////////////////////////////////////////////
    // Try to access all elements on the grid. We don't
    // check whether the values are correct, but rather
    // just whether all of them are set. Since we don't
    // make any assumptions about the types of grids, this
    // test will only be able to indicate on basic memory
    // errors.
    ///////////////////////////////////////////////////////
    real_t sr = 0, sp1 = 0, sp2 = 0;
    for (len_t ir = 0; ir < rg->GetNr(); ir++) {
        sr += rg->GetR(ir) + rg->GetR_f(ir) - rg->GetR_f(ir+1);
        TQS::FVM::MomentumGrid *mg = rg->GetMomentumGrid(ir);

        for (len_t j = 0; j < mg->GetNp2(); j++) {
            sp2 += mg->GetP2(j) + (mg->GetP2_f(j) + mg->GetP2_f(j+1));

            for (len_t i = 0; i < mg->GetNp1(); i++) {
                sp1 += mg->GetP1(i) + (mg->GetP1_f(i) + mg->GetP1_f(i+1));
            }
        }
    }

    return success;
}

/**
 * Run all RadialGrid tests.
 * Returns 'true' if all tests passed. 'false' otherwise.
 */
bool RadialGrid::Run(bool) {
    bool success = true;
    
    // Construct a general grid (where each radius has its
    // own momentum grid)
    if (CheckGeneralGrid())
        this->PrintOK("Successfully constructed a general grid.");
    else {
        success = false;
        this->PrintError("Failed to construct a general copmutational grid.");
    }
    
    return success;
}

