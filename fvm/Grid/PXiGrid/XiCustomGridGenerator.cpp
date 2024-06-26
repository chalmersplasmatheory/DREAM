/**
 * Implementation of a simple uniform mometum (p) grid generator.
 */

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiCustomGridGenerator.hpp"


using namespace DREAM::FVM::PXiGrid;

/**
 * Constructor.
 *
 * xi_f_in: grid points on the pitch flux grid
 *     nxi: Number of points on (cell) grid.
 */
XiCustomGridGenerator::XiCustomGridGenerator(
    const real_t *xi_f_in, const len_t nxi
) : nxi(nxi) {
    this->xif_provided = new real_t[nxi+1];
    for(len_t i=0; i<nxi+1; i++)
        this->xif_provided[i] = xi_f_in[i];
}


/***********************************
 * PUBLIC METHODS                  *
 ***********************************/
/**
 * Re-build the given momentum grid using this
 * grid generator.
 *
 * mg: Momentum grid to re-build.
 *
 * (All other parameters are unused).
 */
bool XiCustomGridGenerator::Rebuild(const real_t, const len_t, MomentumGrid *mg, const RadialGrid*) {
    real_t
        *xi    = new real_t[this->nxi],
        *xi_f  = new real_t[this->nxi+1],
        *dxi   = new real_t[this->nxi],
        *dxi_f = nullptr;

    // Build flux grid
    if(xif_provided==nullptr)
        throw FVMException("XiCustomGridGenerator: Rebuilding grid not yet supported.");
    for (len_t i = 0; i < this->nxi+1; i++)
        xi_f[i] = this->xif_provided[i];
    xif_provided=nullptr;

    for (len_t i = 0; i < this->nxi; i++)
        dxi[i] = xi_f[i+1]-xi_f[i];

    // Build cell grid
    for (len_t i = 0; i < this->nxi; i++)
        xi[i] = 0.5 * (xi_f[i+1] + xi_f[i]);

    if (this->nxi > 1) {
        dxi_f = new real_t[this->nxi-1];
        for (len_t i = 0; i < this->nxi-1; i++)
            dxi_f[i] = xi[i+1] - xi[i];
    }

    mg->InitializeP2("xi", this->nxi, xi, xi_f, dxi, dxi_f);

    initialized = true;
    return true;
}

