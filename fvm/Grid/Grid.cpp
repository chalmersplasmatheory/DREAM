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
 * Get the total number of cells on the radial flux grid,
 * including on the momentum grids at each radius.
 */
const len_t Grid::GetNCells_fr() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells();

    // XXX here we assume that all momentum grids are the same
    // double count last momentum grid (since nr+1)
    N += this->momentumGrids[Nr-1]->GetNCells();

    return N;
}

/**
 * Get the total number of cells on the p1 flux grid.
 */
const len_t Grid::GetNCells_f1() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells_f1();

    return N;
}

/**
 * Get the total number of cells on the p2 flux grid.
 */
const len_t Grid::GetNCells_f2() const {
    const len_t Nr = this->GetNr();
    len_t N = 0;

    for (len_t i = 0; i < Nr; i++)
        N += this->momentumGrids[i]->GetNCells_f2();

    return N;
}

/**
 * Integrate the given vector numerically over the entire
 * phase space (radius+momentum).
 *
 * vec: Vector to integrate numerically.
 */
real_t Grid::Integral(const real_t *vec) const {
    real_t I = 0;

    for (len_t ir = 0, offset=0; ir < this->GetNr(); ir++) {
        real_t VpVol = this->GetVpVol(ir);
        real_t dr    = this->GetRadialGrid()->GetDr(ir);
        I += this->IntegralMomentumAtRadius(ir, vec+offset) * VpVol * dr;

        offset += this->GetMomentumGrid(ir)->GetNCells();
    }

    return I;
}

/**
 * Integrate the given vector numerically over the momentum
 * grids at all radii. The result is stored in the vector
 * 'I', which must be of size 'nr'. The vector 'I' is
 * subsequently returned. If 'I' is a nullptr, a new vector
 * of size 'nr' is allocated and returned.
 *
 * vec: Vector to integrate numerically.
 * I:   Contains integral on return. Must be of size 'nr'. If
 *      'nullptr', this method allocates a new vector of size
 *      'nr', stores the result in it and returns it.
 */
real_t *Grid::IntegralMomentum(const real_t *vec, real_t *I) const {
    if (I == nullptr)
        I = new real_t[this->GetNr()];

    for (len_t ir = 0, offset=0; ir < this->GetNr(); ir++) {
        I[ir] = this->IntegralMomentumAtRadius(ir, vec+offset);
        offset += this->GetMomentumGrid(ir)->GetNCells();
    }

    return I;
}

/**
 * Integrate the given vector numerically over the momentum
 * grid corresponding to the specified radius.
 *
 * ir:  Index of radius to integrate over.
 * vec: Vector to integrate of size np1*np2.
 */
real_t Grid::IntegralMomentumAtRadius(const len_t ir, const real_t *vec) const {
    MomentumGrid *mg = this->GetMomentumGrid(ir);
    const len_t np1 = mg->GetNp1(), np2 = mg->GetNp2();
    const real_t *Vp = this->GetVp(ir);
    const real_t VpVol = this->GetVpVol(ir);

    real_t I = 0;
    for (len_t j = 0; j < np2; j++) {
        real_t dp2 = mg->GetDp2(j);

        for (len_t i = 0; i < np1; i++) {
            real_t dp1 = mg->GetDp1(i);
            len_t idx = j*np1 + i;

            I += vec[idx]*Vp[idx] * dp1*dp2;
        }
    }
    
    return I/VpVol;
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
