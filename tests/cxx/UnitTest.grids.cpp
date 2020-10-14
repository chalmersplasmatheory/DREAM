/**
 * Implementation of various grid generation routines.
 */

#include <string>
#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"
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
/*DREAM::FVM::RadialGrid *UnitTest::InitializeGeneralGridPXi(const len_t nr, const len_t np, const len_t nxi) {
    const real_t B0 = 2;
    const real_t pMin = 0, pMax = 10;

    auto *crgg = new DREAM::FVM::CylindricalRadialGridGenerator(nr, B0);
    auto *rg = new DREAM::FVM::RadialGrid(crgg);

    // Build momentum grid
    for (len_t i = 0; i < nr; i++) {
        auto *pgg = new DREAM::FVM::PXiGrid::PUniformGridGenerator(np+i, pMin, pMax);
        auto *xgg = new DREAM::FVM::PXiGrid::XiUniformGridGenerator(nxi+i);

        auto *mgg = new DREAM::FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
        auto *mg  = new DREAM::FVM::PXiGrid::PXiMomentumGrid(mgg, i, rg);

        rg->SetMomentumGrid(i, mg);
    }

    rg->Rebuild(0);

    return rg;
}*/

/**
 * Initialize a r/p/xi grid with the same momentum grid
 * at all radii.
 *
 * nr:  (optional) Number of radial grid points.
 * np:  (optional) Number of momentum grid points (successively increased by 1 at each radius).
 * nxi: (optional) Number of pitch grid points (sucessively increased by 1 at each radius).
 */
DREAM::FVM::Grid *UnitTest::InitializeGridRCylPXi(
    const len_t nr, const len_t np, const len_t nxi, const real_t B0,
    const real_t pMin, const real_t pMax
) {
    auto *crgg = new DREAM::FVM::CylindricalRadialGridGenerator(nr, B0);
    auto *rg   = new DREAM::FVM::RadialGrid(crgg);

    // Build momentum grid
    auto *pgg = new DREAM::FVM::PXiGrid::PUniformGridGenerator(np, pMin, pMax);
    auto *xgg = new DREAM::FVM::PXiGrid::XiUniformGridGenerator(nxi);

    auto *mgg = new DREAM::FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
    auto *mg  = new DREAM::FVM::PXiGrid::PXiMomentumGrid(mgg, 0, rg);

    auto *grid = new DREAM::FVM::Grid(rg, mg);

    grid->Rebuild(0);

    return grid;
}

/**
 * Initialize a fluid grid (only radial grid).
 *
 * nr: (optional) Number of radial grid points.
 */
DREAM::FVM::Grid *UnitTest::InitializeFluidGrid(const len_t nr, const real_t B0) {
    
    auto *crgg = new DREAM::FVM::CylindricalRadialGridGenerator(nr, B0);
    auto *rg   = new DREAM::FVM::RadialGrid(crgg);

    auto *grid = new DREAM::FVM::Grid(rg, new DREAM::FVM::EmptyMomentumGrid(rg));
    grid->Rebuild(0);

    return grid;
}


DREAM::FVM::Grid *UnitTest::InitializeGridGeneralRPXi(
    const len_t nr, const len_t np, const len_t nxi,
    const len_t ntheta_interp, const len_t nrProfiles, const real_t pMin, const real_t pMax
) {
    real_t r0 = 0.7531;
    real_t ra = 1.4842;
    real_t R0 = 3.39431;
    real_t 
        *rProfiles = new real_t[nrProfiles],
        *Gs = new real_t[nrProfiles], 
        *psi_p0s = new real_t[nrProfiles], 
        *kappas = new real_t[nrProfiles], 
        *deltas = new real_t[nrProfiles], 
        *Deltas = new real_t[nrProfiles];

    for (len_t it = 0; it<nrProfiles; it++){
        rProfiles[it] = r0 + it*(ra-r0)/(nrProfiles-1);
        Gs[it]      = (2.8353 + 2.1828*it*(it-1)/(nrProfiles*nrProfiles))/R0; 
        kappas[it]  = 1.82094 + 1.240021*it/(nrProfiles-1);     
        Deltas[it]  = 0.47183*it/(nrProfiles-1);
        psi_p0s[it] = 1.3913*it/(nrProfiles-1) / R0;        
//        deltas[it]  = 1.9531*it/(nrProfiles-1);
        deltas[it]  = 0.4531*it/(nrProfiles-1);
    }

    struct DREAM::FVM::AnalyticBRadialGridGenerator::shape_profiles shapes = {
        nrProfiles, nrProfiles, nrProfiles, nrProfiles, nrProfiles,
        Gs, rProfiles, psi_p0s, rProfiles, kappas, rProfiles,
        deltas, rProfiles, Deltas, rProfiles
    };

    auto *ABrgg = new DREAM::FVM::AnalyticBRadialGridGenerator(
        nr, r0, ra, R0, ntheta_interp, &shapes
    );

    auto *rg   = new DREAM::FVM::RadialGrid(ABrgg);

    // Construct momentum grid
    auto *pgg = new DREAM::FVM::PXiGrid::PUniformGridGenerator(np, pMin, pMax);
    auto *xgg = new DREAM::FVM::PXiGrid::XiUniformGridGenerator(nxi);

    auto *mgg = new DREAM::FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
    auto *mg  = new DREAM::FVM::PXiGrid::PXiMomentumGrid(mgg, 0, rg);

    // Assemble to full 'Grid' object...
    auto *grid = new DREAM::FVM::Grid(rg, mg);
    grid->Rebuild(0);
    return grid;
}

/**
 * Fluid grid companion to the 'InitializeGridGeneralRPXi()'
 * method above.
 */
DREAM::FVM::Grid *UnitTest::InitializeGridGeneralFluid(
    const len_t nr, const len_t ntheta_interp,
    const len_t nrProfiles
) {
    real_t r0 = 0.7531;
    real_t ra = 1.4842;
    real_t R0 = 3.39431;
    real_t 
        *rProfiles = new real_t[nrProfiles],
        *Gs = new real_t[nrProfiles], 
        *psi_p0s = new real_t[nrProfiles], 
        *kappas = new real_t[nrProfiles], 
        *deltas = new real_t[nrProfiles], 
        *Deltas = new real_t[nrProfiles];

    for (len_t it = 0; it<nrProfiles; it++){
        rProfiles[it] = r0 + it*(ra-r0)/(nrProfiles-1);
        Gs[it]      = (2.8353 + 2.1828*it*(it-1)/(nrProfiles*nrProfiles))/R0; 
        kappas[it]  = 1.82094 + 1.240021*it/(nrProfiles-1);     
        Deltas[it]  = 0.47183*it/(nrProfiles-1);
        psi_p0s[it] = 1.3913*it/(nrProfiles-1);        
        deltas[it]  = 0.4531*it/(nrProfiles-1);
    }

    struct DREAM::FVM::AnalyticBRadialGridGenerator::shape_profiles shapes = {
        nrProfiles, nrProfiles, nrProfiles, nrProfiles, nrProfiles,
        Gs, rProfiles, psi_p0s, rProfiles, kappas, rProfiles,
        deltas, rProfiles, Deltas, rProfiles
    };

    auto *ABrgg = new DREAM::FVM::AnalyticBRadialGridGenerator(
        nr, r0, ra, R0, ntheta_interp, &shapes
    ); 

    auto *rg = new DREAM::FVM::RadialGrid(ABrgg);

    auto *grid = new DREAM::FVM::Grid(rg, new DREAM::FVM::EmptyMomentumGrid(rg));
    grid->Rebuild(0);

    return grid;
}

/**
 * A helper routine for iterating over several different
 * types of grids.
 *
 * igrid: Index of grid to generate.
 */
struct UnitTest::gridcontainer *UnitTest::GetNextGrid(const len_t igrid) {
    string name;
    DREAM::FVM::Grid *rg;

    switch (igrid) {
        case 0:
            //rg = InitializeGeneralGridPXi();
            rg = InitializeGridRCylPXi();
            name = "Cylindrical r, single p/xi grid";
            break;

        default: return nullptr;
    }

    struct UnitTest::gridcontainer *gc = new UnitTest::gridcontainer;
    gc->name = name;
    gc->grid = rg;

    return gc;
}

