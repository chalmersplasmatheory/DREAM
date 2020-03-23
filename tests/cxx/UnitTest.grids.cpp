/**
 * Implementation of various grid generation routines.
 */

#include <string>
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "UnitTest.hpp"


using namespace std;
using namespace DREAMTESTS;

/**
 * Initialize a default grid.
 *
 * nr:  (optional) Number of radial grid points.
 * np:  (optional) Number of momentum grid points (successively increased by 1 at each radius).
 * nxi: (optional) Number of pitch grid points (sucessively increased by 1 at each radius).
 */
DREAM::FVM::RadialGrid *UnitTest::InitializeGeneralGridPXi(const len_t nr, const len_t np, const len_t nxi) {
    const real_t B0 = 2;
    const real_t pMin = 0, pMax = 10;

    auto *crgg = new DREAM::FVM::CylindricalRadialGridGenerator(nr, B0);
    auto *rg = new DREAM::FVM::RadialGrid(crgg);

    // Build momentum grid
    for (len_t i = 0; i < nr; i++) {
        auto *pgg = new DREAM::FVM::PXiGrid::PUniformGridGenerator(np+i, pMin, pMax);
        auto *xgg = new DREAM::FVM::PXiGrid::XiUniformGridGenerator(nxi+i);

        auto *mgg = new DREAM::FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
        auto *mg  = new DREAM::FVM::MomentumGrid(mgg, i, rg);

        rg->SetMomentumGrid(i, mg);
    }

    rg->Rebuild(0);

    return rg;
}

/**
 * A helper routine for iterating over several different
 * types of grids.
 *
 * igrid: Index of grid to generate.
 */
struct UnitTest::gridcontainer *UnitTest::GetNextGrid(const len_t igrid) {
    string name;
    DREAM::FVM::RadialGrid *rg;

    switch (igrid) {
        case 0:
            rg = InitializeGeneralGridPXi();
            name = "General p/xi grid";
            break;

        default: return nullptr;
    }

    struct UnitTest::gridcontainer *gc = new UnitTest::gridcontainer;
    gc->name = name;
    gc->grid = rg;

    return gc;
}

