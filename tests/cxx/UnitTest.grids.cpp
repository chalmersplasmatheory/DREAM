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
#include "FVM/Grid/NumericBRadialGridGenerator.hpp"
#include "FVM/Grid/Stellarator/NumericStellaratorRadialGridGenerator.hpp"
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

/**
 * Initialize a fluid grid (only radial grid).
 *
 * nr: (optional) Number of radial grid points.
 */
DREAM::FVM::Grid *UnitTest::InitializeNumericFluidGrid(const real_t *r_f, const len_t nr, const std::string& mf) {

    auto *NBrgg = new DREAM::FVM::NumericBRadialGridGenerator(r_f, nr, mf);

    auto *rg   = new DREAM::FVM::RadialGrid(NBrgg);

    auto *grid = new DREAM::FVM::Grid(rg, new DREAM::FVM::EmptyMomentumGrid(rg));
    grid->Rebuild(0);

    return grid;
}

/**
 * Initialize a fluid grid (only radial grid).
 *
 * nr: (optional) Number of radial grid points.
 */
DREAM::FVM::Grid *UnitTest::InitializeNumericStellaratorFluidGrid(const real_t *r_f, const len_t nr, SFile_HDF5 *s) {
    real_t R0 = s->GetDoubles1D("R0", nullptr)[0];
    len_t nfp = s->GetInt64_1D("nfp", nullptr)[0];
    len_t ntheta_interp = s->GetInt64_1D("ntheta", nullptr)[0];
    len_t nphi_interp   = s->GetInt64_1D("nphi", nullptr)[0];
    real_t b  = s->GetDoubles1D("wall_radius", nullptr)[0];

    DREAM::FVM::NumericStellaratorRadialGridGenerator::eq_data *eqdata =
        new DREAM::FVM::NumericStellaratorRadialGridGenerator::eq_data;
    
    sfilesize_t dims; 
    eqdata->rho             = s->GetDoubles1D("rho", &dims);
    eqdata->nrho            = dims;
    eqdata->theta           = s->GetDoubles1D("theta", &dims);
    eqdata->ntheta          = dims;
    eqdata->phi             = s->GetDoubles1D("phi", &dims);
    eqdata->nphi            = dims;
    eqdata->dataR           = s->GetDoubles1D("R", nullptr);
    eqdata->dataZ           = s->GetDoubles1D("Z", nullptr);
    eqdata->dataG           = s->GetDoubles1D("G", nullptr);
    eqdata->dataI           = s->GetDoubles1D("I", nullptr);
    eqdata->dataiota        = s->GetDoubles1D("iota", nullptr);
    eqdata->datapsi         = s->GetDoubles1D("psi_T", nullptr);
    eqdata->dataB           = s->GetDoubles1D("B", nullptr);
    eqdata->dataBdotGradphi = s->GetDoubles1D("BdotGradPhi", nullptr);
    eqdata->dataJacobian    = s->GetDoubles1D("Jacobian", nullptr);
    eqdata->datagtt         = s->GetDoubles1D("g_tt", nullptr);
    eqdata->datagtp         = s->GetDoubles1D("g_tp", nullptr);
    eqdata->datalambdat     = s->GetDoubles1D("lambda_t", nullptr);
    eqdata->datalambdap     = s->GetDoubles1D("lambda_p", nullptr);

    auto *NBrgg = new DREAM::FVM::NumericStellaratorRadialGridGenerator(r_f, nr, b, R0, nfp, eqdata, ntheta_interp, nphi_interp);

    auto *rg   = new DREAM::FVM::RadialGridStellarator(NBrgg);

    auto *grid = new DREAM::FVM::Grid(rg, new DREAM::FVM::EmptyMomentumGrid(rg));
    grid->Rebuild(0);

    return grid;
}


DREAM::FVM::Grid *UnitTest::InitializeGridGeneralRPXi(
    const len_t nr, const len_t np, const len_t nxi,
    const len_t ntheta_interp, const len_t nrProfiles, const real_t pMin, const real_t pMax,
    DREAM::FVM::FluxSurfaceAverager::quadrature_method q_method_passing,
    DREAM::FVM::FluxSurfaceAverager::quadrature_method q_method_trapped
) {
    real_t r0 = 0;
    real_t ra = 2;
    real_t R0 = 4;

    real_t 
        DeltaMax = 0.6,
        deltaMax = 0.2,
        GMin = 4,
        GMax = 4.5,
        kappaMin = 1.4,
        kappaMax = 1.9,
        DeltaPsi = M_PI;
        

    real_t 
        *rProfiles = new real_t[nrProfiles],
        *Gs = new real_t[nrProfiles], 
        *psi_p0s = new real_t[nrProfiles], 
        *kappas = new real_t[nrProfiles], 
        *deltas = new real_t[nrProfiles], 
        *Deltas = new real_t[nrProfiles];

    #define ProfileAtIt(Y,YMIN,YMAX) Y[it] = YMIN + (YMAX-YMIN)*it/(nrProfiles-1)
    for (len_t it = 0; it<nrProfiles; it++){
        ProfileAtIt(rProfiles, r0, ra);
        ProfileAtIt(Gs, GMin, GMax);
        ProfileAtIt(kappas, kappaMin, kappaMax);
        ProfileAtIt(Deltas, 0, DeltaMax);
        ProfileAtIt(deltas, 0, deltaMax);
        psi_p0s[it] = DeltaPsi*it*it/((nrProfiles-1)*(nrProfiles-1));
    }
    #undef ProfileAtIt

    struct DREAM::FVM::AnalyticBRadialGridGenerator::shape_profiles shapes = {
        nrProfiles, nrProfiles, nrProfiles, nrProfiles, nrProfiles,
        Gs, rProfiles, psi_p0s, rProfiles, kappas, rProfiles,
        deltas, rProfiles, Deltas, rProfiles
    };

    auto *ABrgg = new DREAM::FVM::AnalyticBRadialGridGenerator(
        nr, r0, ra, R0, ntheta_interp, &shapes
    );

    auto *rg   = new DREAM::FVM::RadialGrid(ABrgg, 0, DREAM::FVM::FluxSurfaceAverager::INTERP_STEFFEN, q_method_passing);

    // Construct momentum grid
    auto *pgg = new DREAM::FVM::PXiGrid::PUniformGridGenerator(np, pMin, pMax);
    auto *xgg = new DREAM::FVM::PXiGrid::XiUniformGridGenerator(nxi);

    auto *mgg = new DREAM::FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
    auto *mg  = new DREAM::FVM::PXiGrid::PXiMomentumGrid(mgg, 0, rg);

    // Assemble to full 'Grid' object...
    auto *grid = new DREAM::FVM::Grid(rg, mg, 0, q_method_trapped);
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
) {    real_t r0 = 0;
    real_t ra = 2;
    real_t R0 = 4;

    real_t 
        DeltaMax = 0.6,
        deltaMax = 0.2,
        GMin = 4,
        GMax = 4.5,
        kappaMin = 1.4,
        kappaMax = 1.9,
        DeltaPsi = M_PI;
        

    real_t 
        *rProfiles = new real_t[nrProfiles],
        *Gs = new real_t[nrProfiles], 
        *psi_p0s = new real_t[nrProfiles], 
        *kappas = new real_t[nrProfiles], 
        *deltas = new real_t[nrProfiles], 
        *Deltas = new real_t[nrProfiles];

    #define ProfileAtIt(Y,YMIN,YMAX) Y[it] = YMIN + (YMAX-YMIN)*it/(nrProfiles-1)
    for (len_t it = 0; it<nrProfiles; it++){
        ProfileAtIt(rProfiles, r0, ra);
        ProfileAtIt(Gs, GMin, GMax);
        ProfileAtIt(kappas, kappaMin, kappaMax);
        ProfileAtIt(Deltas, 0, DeltaMax);
        ProfileAtIt(deltas, 0, deltaMax);
        psi_p0s[it] = DeltaPsi*it*it/((nrProfiles-1)*(nrProfiles-1));
    }
    #undef ProfileAtIt
    
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
        case 1:
            rg = InitializeGridGeneralRPXi();
            name = "General r, single p/xi grid";
            break;
        default: return nullptr;
    }

    struct UnitTest::gridcontainer *gc = new UnitTest::gridcontainer;
    gc->name = name;
    gc->grid = rg;

    return gc;
}

