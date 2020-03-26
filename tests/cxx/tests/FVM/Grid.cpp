/**
 * Test for the 'Grid' class in the FVM library.
 */

#include "UnitTest.hpp"
#include "Grid.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Check whether a general radial grid can be appropriately
 * constructed.
 */
bool Grid::CheckGridRCylPXi() {
    bool success = true;

    DREAM::FVM::Grid *grid = this->InitializeGridRCylPXi();

    ///////////////////////////////////////////////////////
    // Try to access all elements on the grid. We don't
    // check whether the values are correct, but rather
    // just whether all of them are set. Since we don't
    // make any assumptions about the types of grids, this
    // test will only be able to indicate on basic memory
    // errors.
    ///////////////////////////////////////////////////////
    real_t sr = 0, sp1 = 0, sp2 = 0;
    for (len_t ir = 0; ir < grid->GetNr(); ir++) {
        sr += grid->GetRadialGrid()->GetR(ir)
            + grid->GetRadialGrid()->GetR_f(ir)
            - grid->GetRadialGrid()->GetR_f(ir+1);

        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);

        for (len_t j = 0; j < mg->GetNp2(); j++) {
            sp2 += mg->GetP2(j) + (mg->GetP2_f(j) + mg->GetP2_f(j+1));

            for (len_t i = 0; i < mg->GetNp1(); i++) {
                sp1 += mg->GetP1(i) + (mg->GetP1_f(i) + mg->GetP1_f(i+1));
            }
        }
    }

    delete grid;

    return success;
}

/**
 * Run all Grid tests.
 * Returns 'true' if all tests passed. 'false' otherwise.
 */
bool Grid::Run(bool) {
    bool success = true;
    
    // Construct a general grid (where each radius has its
    // own momentum grid)
    if (CheckGridRCylPXi())
        this->PrintOK("Successfully constructed a cylindrical r-grid, with p/xi momentum grid.");
    else {
        success = false;
        this->PrintError("Failed to construct a cylindrical r-grid with p/xi momentum grid.");
    }

    // TODO
    // - Multiple radii, one momentum grid
    // - One radius, full momentum grid
    
    return success;
}

