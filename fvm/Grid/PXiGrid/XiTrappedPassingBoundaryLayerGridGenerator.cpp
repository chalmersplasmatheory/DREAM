/**
 * This grid automatically sets up a xi grid based on the location
 * of the trapped-passing boundaries at each radius.
 */

#include <algorithm>
#include <vector>
#include "FVM/Grid/PXiGrid/XiTrappedPassingBoundaryLayerGridGenerator.hpp"


using namespace DREAM::FVM::PXiGrid;
using namespace std;

/**
 * Constructor.
 *
 * dxiMax:             Maximum allowed grid spacing dxi: if this spacing is
 *                     exceeded after placing points on the trapped-passing
 *                     boundary, will fill in the gap by adding one or more grid
 *                     points (uniformly) in the gaps.
 * nxiPass:            Number of grid cells in pitch to be contained between the
 *                     largest xi0Trapped and xi0=+/-1.
 * nxiTrap:            Number of grid cells in pitch to be contained between the
 *                     (in absolute value) minimum xi0Trapped and xi0.
 * boundaryLayerWidth: Width in pitch of the grid cell containing each
 *                     trapped-passing boundary, typically << 1.
 */
XiTrappedPassingBoundaryLayerGridGenerator::XiTrappedPassingBoundaryLayerGridGenerator(
    const real_t dxiMax, const len_t nxiPass, const len_t nxiTrap,
    const real_t boundaryLayerWidth
) : dxiMax(dxiMax), nxiPass(nxiPass), nxiTrap(nxiTrap),
    boundaryLayerWidth(boundaryLayerWidth) {
}

/**
 * Re-build the given momentum grid using this grid generator.
 *
 * mg: Momentum grid to re-build.
 * rg: Radial grid on which this momentum grid lives.
 */
bool XiTrappedPassingBoundaryLayerGridGenerator::Rebuild(
    const real_t, const len_t, MomentumGrid *mg, const RadialGrid *rg
) {
    const len_t nr = rg->GetNr();
    const real_t *xi0Trapped = rg->GetXi0TrappedBoundary();
    vector<real_t> xi_f;

    // Locate max and min boundaries...
    real_t xi0Trapped_max=0, xi0Trapped_min=1;
    for (len_t ir = 0; ir < nr; ir++) {
        if (xi0Trapped[ir] > xi0Trapped_max)
            xi0Trapped_max = xi0Trapped[ir];
        if (xi0Trapped[ir] < xi0Trapped_min)
            xi0Trapped_min = xi0Trapped[ir];
    }

    // Insert passing region points [1, +xi0Trapped) and (-xi0Trapped, -1]...
    for (len_t i = 0; i < nxiPass; i++) {
        real_t x = 1 - i*(1-xi0Trapped_max) / nxiPass;
        xi_f.push_back(+x);
        xi_f.push_back(-x);
    }

    // Insert lowest part of trapped region (if any)
    if (nxiTrap > 1 && xi0Trapped_min > 0) {
        for (len_t i = 1; i < nxiTrap; i++) {
            real_t x = i * xi0Trapped_min / nxiTrap;
            xi_f.push_back(x);
            xi_f.push_back(-x);
        }
    }

    // Add xi=0
    xi_f.push_back(0.0);

    // Add cells around each trapped-passing boundary...
    for (len_t ir = 0; ir < nr; ir++) {
        if (xi0Trapped[ir] > 0) {
            real_t xiAdd1 = xi0Trapped[ir] + 0.5*this->boundaryLayerWidth;
            real_t xiAdd2 = xi0Trapped[ir] - 0.5*this->boundaryLayerWidth;

            xi_f.push_back(-xiAdd1);
            xi_f.push_back(xiAdd1);
            xi_f.push_back(-xiAdd2);
            xi_f.push_back(xiAdd2);
        }
    }

    std::sort(xi_f.begin(), xi_f.end());

    // Add additional points if the resulting spacing
    // exceeds the desired dxiMax...
    for (vector<real_t>::iterator it = xi_f.begin()+1; it != xi_f.end(); it++) {
        real_t dxi = *it - *(it-1);
        int_t N = (int_t)std::ceil(dxi / this->dxiMax);

        if (N > 1) {
            real_t x1 = *(it-1), x2 = *it;
            it -= 1;
            for (int_t i = 1; i < N+1; i++)
                it = xi_f.insert(it+1, x1 + i*(x2-x1)/(N+1));

            it++;
        }
    }

    // Set arrays...
    this->nxi = xi_f.size()-1;
    real_t
        *_xi   = new real_t[this->nxi],
        *_xi_f = new real_t[this->nxi+1],
        *dxi   = new real_t[this->nxi],
        *dxi_f = nullptr;

    for (len_t i = 0; i < this->nxi+1; i++)
        _xi_f[i] = xi_f[i];
    for (len_t i = 0; i < this->nxi; i++)
        dxi[i] = _xi_f[i+1] - _xi_f[i];
    for (len_t i = 0; i < this->nxi; i++)
        _xi[i] = 0.5*(_xi_f[i+1]+_xi_f[i]);

    if (this->nxi > 1) {
        dxi_f = new real_t[this->nxi-1];
        for (len_t i = 0; i < this->nxi-1; i++)
            dxi_f[i] = _xi[i+1] - _xi[i];
    }

    mg->InitializeP2("xi", this->nxi, _xi, _xi_f, dxi, dxi_f);

    initialized = true;
    return true;
}

