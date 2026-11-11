/**
 * Test for the 'Grid' class in the FVM library.
 */

#include "UnitTest.hpp"
#include "Grid.hpp"

#include <cmath>
#include <memory>
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiCustomGridGenerator.hpp"


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
 * Compare the metric in a cell that straddles zero with the direct pitch
 * integral of the represented trapped orbit.
 */
bool Grid::CheckZeroCrossingPitchCell() {
    using namespace DREAM::FVM;
    static const real_t r[] = {0, 1}, one[] = {1, 1}, zero[] = {0, 0};
    auto *profiles = new AnalyticBRadialGridGenerator::shape_profiles{
        2, 2, 2, 2, 2, one, r, zero, r, one, r, zero, r, zero, r
    };
    auto radial = std::make_unique<RadialGrid>(
        new AnalyticBRadialGridGenerator(1, 0.29, 0.31, 1, 128, profiles)
    );
    const real_t faces[] = {-1, -0.01, 0.02, 1};
    auto *momentum = new PXiGrid::MomentumGridGenerator(
        new PXiGrid::PUniformGridGenerator(1, 0, 2),
        new PXiGrid::XiCustomGridGenerator(faces, 3)
    );
    DREAM::FVM::Grid grid(radial.get(), new PXiGrid::PXiMomentumGrid(momentum, 0, radial.get()));
    grid.Rebuild(0);

    const real_t reference = radial->GetFluxSurfaceAverager()->EvaluateCellAveragedBounceIntegralOverP2(
        0, faces[1], faces[2], FLUXGRIDTYPE_DISTRIBUTION,
        RadialGrid::BA_FUNC_UNITY, nullptr, RadialGrid::BA_PARAM_UNITY
    );
    const real_t metric = grid.GetVpOverP2AtZero(0)[1];
    // The grid uses a linear approximation in pitch for this narrow cell.
    if (!(reference > 0 && std::abs(metric - reference) < 5e-3 * reference)) {
        this->PrintError("Zero-crossing pitch cell metric %.8g, direct average %.8g.", metric, reference);
        return false;
    }
    return true;
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

    if (CheckZeroCrossingPitchCell())
        this->PrintOK("The zero-crossing pitch cell metric agrees with the direct average.");
    else
        success = false;

    // TODO
    // - Multiple radii, one momentum grid
    // - One radius, full momentum grid
    
    return success;
}

