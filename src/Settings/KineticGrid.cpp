/**
 * Routines for constructing kinetic (i.e. hot-tail and runaway)
 * grids.
 */

#include <string>
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PBiUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PCustomGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiCustomGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiBiUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiBiUniformThetaGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiUniformThetaGridGenerator.hpp"


using namespace DREAM;
using namespace std;

// Module names (to cause compile-time errors if
// misspelled, instead of run-time errors)
#define HOTTAILGRID "hottailgrid"
#define RUNAWAYGRID "runawaygrid"


/**
 * Define options common to all kinetic grids. This method
 * is called by 'DefineOptions_HotTailGrid()' and
 * 'DefineOptions_RunawayGrid()'.
 *
 * mod: Name of the module owning the setting.
 * s:   Settings object to define settings in.
 */
void SimulationGenerator::DefineOptions_KineticGrid(const string& mod, Settings *s) {
    s->DefineSetting(mod + "/enabled", "Indicates whether this momentum grid is used in the simulation", (bool)false, true);
    s->DefineSetting(mod + "/type", "Momentum grid type", (int_t)OptionConstants::MOMENTUMGRID_TYPE_PXI);

    // p/xi grid
    s->DefineSetting(mod + "/np", "Number of distribution grid points in p", (int_t)1);
    s->DefineSetting(mod + "/nxi", "Number of distribution grid points in xi", (int_t)1);
    s->DefineSetting(mod + "/pmax", "Maximum momentum on the (flux) grid", (real_t)0.0);
    s->DefineSetting(mod + "/pgrid", "Type of momentum grid to generate", (int_t)OptionConstants::PXIGRID_PTYPE_UNIFORM);
    s->DefineSetting(mod + "/xigrid", "Type of pitch grid to generate", (int_t)OptionConstants::PXIGRID_XITYPE_UNIFORM);

    // nonuniform p grid
    s->DefineSetting(mod + "/npsep", "Number of distribution grid points for pmin<p<psep", (int_t)0);
    s->DefineSetting(mod + "/npsep_frac", "Fraction of distribution grid points for pmin<p<psep", (real_t)0.0);
    s->DefineSetting(mod + "/psep", "Separating momentum on the biuniform (flux) grid", (real_t)0.0);
    
    // nonuniform xi grid
    s->DefineSetting(mod + "/nxisep", "Number of distribution grid points for xisep<xi<1", (int_t)0);
    s->DefineSetting(mod + "/nxisep_frac", "Fraction of distribution grid points for xisep<xi<1", (real_t)0.0);
    s->DefineSetting(mod + "/xisep", "Separating pitch on the biuniform (flux) grid", (real_t)-1.0);

    // custom grids
    s->DefineSetting(mod + "/p_f",  "Grid points of the momentum flux grid", 0, (real_t*) nullptr);
    s->DefineSetting(mod + "/xi_f", "Grid points of the pitch flux grid", 0, (real_t*) nullptr);
   
}

/**
 * Define options applicable to the hot-tail grid.
 *
 * s:   Settings object to define settings in.
 */
void SimulationGenerator::DefineOptions_HotTailGrid(Settings *s) {
    DefineOptions_KineticGrid(HOTTAILGRID, s);
}

/**
 * Define options applicable to the runaway grid.
 *
 * s:   Settings object to define settings in.
 */
void SimulationGenerator::DefineOptions_RunawayGrid(Settings *s) {
    DefineOptions_KineticGrid(RUNAWAYGRID, s);
}


/*******************************
 * GRID CONSTRUCTION           *
 *******************************/
/**
 * Construct the hot-tail grid according to the
 * given specification 's'.
 *
 * s:     Settings object specifying how to construct the
 *        hot-tail grid.
 * rgrid: Radial grid to use for defining the hot-tail grid.
 * type:  On return, contains the exact type of the constructed
 *        momentum grid.
 */
FVM::Grid *SimulationGenerator::ConstructHotTailGrid(
    Settings *s, FVM::RadialGrid *rgrid,
    enum OptionConstants::momentumgrid_type *type
) {
    bool enabled = s->GetBool(HOTTAILGRID "/enabled");

    if (!enabled)
        return nullptr;

    *type = (enum OptionConstants::momentumgrid_type)s->GetInteger(HOTTAILGRID "/type");

    FVM::MomentumGrid *mg;
    switch (*type) {
        case OptionConstants::MOMENTUMGRID_TYPE_PXI:
            mg = Construct_PXiGrid(s, HOTTAILGRID, 0.0, rgrid);
            break;

        // XXX WARNING: The runaway grid assumes that the first coordinate
        // on this grid is 'p'!

        default:
            throw SettingsException(
                "Unrecognized momentum grid type specified to hot-tail grid: " INT_T_PRINTF_FMT ".",
                *type
            );
    }

    return new FVM::Grid(rgrid, mg);
}

/**
 * Construct the runaway grid according to the
 * given specification 's'.
 *
 * s:           Settings object specifying how to construct the
 *              runaway grid.
 * rgrid:       Radial grid to use for defining the runaway grid.
 * hottailGrid: Hot-tail grid to use for defining the runaway grid
 *              (can be 'nullptr').
 * type:        On return, contains the exact type of the constructed
 *              momentum grid.
 */
FVM::Grid *SimulationGenerator::ConstructRunawayGrid(
    Settings *s, FVM::RadialGrid *rgrid, FVM::Grid *hottailGrid,
    enum OptionConstants::momentumgrid_type *type
) {
    bool enabled = s->GetBool(RUNAWAYGRID "/enabled");

    if (!enabled)
        return nullptr;

    *type = (enum OptionConstants::momentumgrid_type)s->GetInteger(RUNAWAYGRID "/type");

    FVM::MomentumGrid *mg;
    real_t pmin;
    switch (*type) {
        case OptionConstants::MOMENTUMGRID_TYPE_PXI:
            if (hottailGrid == nullptr)
                pmin = 0;
            else
                pmin = hottailGrid->GetMomentumGrid(0)->GetP1(
                    hottailGrid->GetMomentumGrid(0)->GetNp1()-1
                );
            mg = Construct_PXiGrid(s, RUNAWAYGRID, pmin, rgrid);
            break;

        default:
            throw SettingsException(
                "Unrecognized momentum grid type specified to runaway grid: " INT_T_PRINTF_FMT ".",
                *type
            );
    }

    return new FVM::Grid(rgrid, mg);
}

/**
 * Construct a p/xi momentum grid.
 *
 * s:    Settings object specifying how to construct the grid.
 * mod:  Name of the module to load settings from.
 * pmin: Minimum momentum value on (flux) grid.
 */
FVM::PXiGrid::PXiMomentumGrid *SimulationGenerator::Construct_PXiGrid(
    Settings *s, const string& mod, const real_t pmin,
    FVM::RadialGrid *rgrid
) {
    int_t  np   = s->GetInteger(mod + "/np");
    int_t  nxi  = s->GetInteger(mod + "/nxi");
    real_t pmax = s->GetReal(mod + "/pmax");

    enum OptionConstants::pxigrid_ptype pgrid   = (enum OptionConstants::pxigrid_ptype)s->GetInteger(mod+"/pgrid");
    enum OptionConstants::pxigrid_xitype xigrid = (enum OptionConstants::pxigrid_xitype)s->GetInteger(mod+"/xigrid");

    // Verify grid limits
    if (pmax <= pmin)
        throw SettingsException(
            "%s: PMAX must be strictly greater than PMIN.",
            mod.c_str()
        );
    
    FVM::PXiGrid::PGridGenerator *pgg;
    FVM::PXiGrid::XiGridGenerator *xgg;

    // Construct P grid generator
    switch (pgrid) {
        case OptionConstants::PXIGRID_PTYPE_UNIFORM:
            pgg = new FVM::PXiGrid::PUniformGridGenerator(np, pmin, pmax);
            break;

        case OptionConstants::PXIGRID_PTYPE_BIUNIFORM: {
            real_t psep = s->GetReal(mod + "/psep");
            len_t npSep = s->GetInteger(mod + "/npsep");
            real_t npSep_frac = s->GetReal(mod + "/npsep_frac");

            if (npSep > 0)
                pgg = new FVM::PXiGrid::PBiUniformGridGenerator(np, (len_t)npSep, pmin, psep, pmax);
            else if (npSep_frac > 0)
                pgg = new FVM::PXiGrid::PBiUniformGridGenerator(np, (real_t)npSep_frac, pmin, psep, pmax);
            else
                throw SettingsException("%s: Neither 'npsep' nor 'npsep_frac' have been specified.", mod.c_str());
        } break;

        case OptionConstants::PXIGRID_PTYPE_CUSTOM:{
            len_t len_pf;
            const real_t *p_f = s->GetRealArray(mod + "/p_f", 1, &len_pf);

            if (p_f[0] != pmin)
                throw SettingsException("%s: The first point on the custom momentum grid must be %f.", mod.c_str(), pmin);

            pgg = new FVM::PXiGrid::PCustomGridGenerator(p_f, len_pf-1);
        } break;

        default:
            throw SettingsException(
                "%s: Unrecognized P grid type specified: %d.",
                mod.c_str(), pgrid
            );
    }

    // Construct XI grid generator
    switch (xigrid) {
        case OptionConstants::PXIGRID_XITYPE_UNIFORM:
            xgg = new FVM::PXiGrid::XiUniformGridGenerator(nxi);
            break;

        case OptionConstants::PXIGRID_XITYPE_BIUNIFORM:{ 
            real_t xisep = s->GetReal(mod + "/xisep");
            len_t nxiSep = s->GetInteger(mod + "/nxisep");
            real_t nxiSep_frac = s->GetReal(mod + "/nxisep_frac");

            if (nxiSep > 0)
                xgg = new FVM::PXiGrid::XiBiUniformGridGenerator(nxi, (len_t)nxiSep, xisep);
            else if (nxiSep_frac > 0)
                xgg = new FVM::PXiGrid::XiBiUniformGridGenerator(nxi, (real_t)nxiSep_frac, xisep);
            else
                throw SettingsException("%s: Neither 'nxisep' nor 'nxisep_frac' have been specified.", mod.c_str());
        } break;
        
        case OptionConstants::PXIGRID_XITYPE_UNIFORM_THETA:
            xgg = new FVM::PXiGrid::XiUniformThetaGridGenerator(nxi);
            break;
        	
        case OptionConstants::PXIGRID_XITYPE_BIUNIFORM_THETA:{ 
            real_t xisep = s->GetReal(mod + "/xisep");
            len_t nxiSep = s->GetInteger(mod + "/nxisep");
            real_t nxiSep_frac = s->GetReal(mod + "/nxisep_frac");

            if (nxiSep > 0)
                xgg = new FVM::PXiGrid::XiBiUniformThetaGridGenerator(nxi, (len_t)nxiSep, xisep);
            else if (nxiSep_frac > 0)
                xgg = new FVM::PXiGrid::XiBiUniformThetaGridGenerator(nxi, (real_t)nxiSep_frac, xisep);
            else
                throw SettingsException("%s: Neither 'nxisep' nor 'nxisep_frac' have been specified.", mod.c_str());
        } break;
	
        case OptionConstants::PXIGRID_XITYPE_CUSTOM:{
            len_t len_xif;
            const real_t *xi_f = s->GetRealArray(mod + "/xi_f", 1, &len_xif);
            xgg = new FVM::PXiGrid::XiCustomGridGenerator(xi_f, len_xif-1);
        } break;
        default:
            throw SettingsException(
                "%s: Unrecognized XI grid type specified: %d.",
                mod.c_str(), xigrid
            );
    }

    auto *pxmgg = new FVM::PXiGrid::MomentumGridGenerator(pgg, xgg);
    return new FVM::PXiGrid::PXiMomentumGrid(pxmgg, 0, rgrid);
}

