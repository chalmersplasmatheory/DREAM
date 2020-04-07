/**
 * Implementation of an overarching 'Grid' object.
 */

#include <algorithm>
#include <vector>
#include "FVM/Grid/Grid.hpp"


using namespace std;
using namespace DREAM::FVM;

/**
 * Constructor.
 */
Grid::Grid(RadialGrid *rg, MomentumGrid *mg, const real_t /*t0*/) {
    this->rgrid = rg;
    this->momentumGrids = new MomentumGrid*[rgrid->GetNr()];

    for (len_t i = 0; i < rgrid->GetNr(); i++)
        this->momentumGrids[i] = mg;
}

/**
 * Destructor.
 */
Grid::~Grid() {
    const len_t nr = this->GetNr();

    // Destroy momentum grids
    //   Since several, or even all, radii may share
    //   a single momentum grid, we should be careful
    //   not to try to double-free any momentum grid.
    vector<MomentumGrid*> deletedPtrs(nr);
    for (len_t i = 0; i < nr; i++) {
        MomentumGrid *p = this->momentumGrids[i];

        // Has the MomentumGrid been deleted already?
        if (find(deletedPtrs.begin(), deletedPtrs.end(), p) != deletedPtrs.end())
            continue;

        deletedPtrs.push_back(p);
        delete p;
    }

    delete [] this->momentumGrids;
    delete this->rgrid;
}

/*****************************
 * PUBLIC METHODS            *
 *****************************/
/**
 * Get the total number of cells on this grid,
 * including on the momentum grids at each radius.
 */
const len_t Grid::GetNCells() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells();

    return N;
}

/**
 * Rebuilds any non-static (i.e. time dependent) grids
 * used. This can be used if, for example, a dynamically
 * evolving magnetic equilibrium is used, or if some
 * grids are adaptive.
 *
 * t: Time to which re-build the grids for.
 */
bool Grid::Rebuild(const real_t t) {
    bool rgridUpdated = false, updated = false;

    // Re-build radial grid
    if (this->rgrid->NeedsRebuild(t))
        rgridUpdated = this->rgrid->Rebuild(t);

    updated = rgridUpdated;

    // Re-build momentum grids
    const len_t Nr = GetNr();
    for (len_t i = 0; i < Nr; i++) {
        if (this->momentumGrids[i]->NeedsRebuild(t, rgridUpdated))
            updated |= this->momentumGrids[i]->Rebuild(t, i, this->rgrid);
    }

    // Re-build jacobians
    if (updated)
        this->RebuildJacobians();

    return updated;
}
