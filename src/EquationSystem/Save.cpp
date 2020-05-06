/**
 * Routines for saving things from the EquationSystem object.
 */

#include <string>
#include <softlib/SFile.h>
#include "DREAM/EquationSystem.hpp"


using namespace DREAM;
using namespace std;


/**
 * Save time, radial and momentum grids as raw vectors
 * to the given SFile object.
 *
 * sf:   SFile object to write grids to.
 * path: Path to group which the data should be written to.
 *       Note that the group must exist in the file.
 */
void EquationSystem::SaveGrids(SFile *sf, const string& path) {
    string group;

    if (path.back() == '/')
        group = path;
    else
        group = path + "/";

    // Time grid
    const real_t *t = this->times.data();
    sf->WriteList(group + "t", t, this->times.size());

    // Radial grid
    const real_t *r = this->fluidGrid->GetRadialGrid()->GetR();
    sf->WriteList(group + "r", r, this->fluidGrid->GetNr());

    // Hot-tail grid
    if (this->hottailGrid != nullptr) {
        sf->CreateStruct(group + "hottail");
        SaveMomentumGrid(sf, group + "hottail/", this->hottailGrid, this->hottailGrid_type);
    }

    // Runaway grid
    if (this->runawayGrid != nullptr) {
        sf->CreateStruct(group + "hottail");
        SaveMomentumGrid(sf, group + "hottail/", this->runawayGrid, this->runawayGrid_type);
    }
}

/**
 * XXX Here we assume that all momentum grids are the same
 * Save a momentum grid to the given SFile object.
 *
 * sf:       SFile object to write grids to.
 * gridpath: Full path to grid data in output file.
 * g:        Grid to save.
 */
void EquationSystem::SaveMomentumGrid(
    SFile *sf, const string& gridname, FVM::Grid *g,
    enum OptionConstants::momentumgrid_type tp
) {
    const real_t *p1 = g->GetMomentumGrid(0)->GetP1();
    const real_t *p2 = g->GetMomentumGrid(0)->GetP2();
    const len_t np1  = g->GetMomentumGrid(0)->GetNp1();
    const len_t np2  = g->GetMomentumGrid(0)->GetNp2();

    // Write grid type
    sf->WriteInt32List(gridname + "type", (int32_t*)&tp, 1);

    sf->WriteList(gridname + "p1", p1, np1);
    sf->WriteList(gridname + "p2", p2, np2);
}

