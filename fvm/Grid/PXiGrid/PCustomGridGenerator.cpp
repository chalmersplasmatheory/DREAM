/**
 * Implementation of a simple uniform mometum (p) grid generator.
 */

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PCustomGridGenerator.hpp"


using namespace DREAM::FVM::PXiGrid;

/**
 * Constructor.
 *
 * p_f_in: grid points on the momentum flux grid
 *     np: Number of points on (cell) grid.
 */
PCustomGridGenerator::PCustomGridGenerator(
    const real_t *p_f_in, const len_t np
) : np(np) {
    this->pf_provided = new real_t[np+1];
    for(len_t i=0; i<np+1; i++)
        this->pf_provided[i] = p_f_in[i];
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
bool PCustomGridGenerator::Rebuild(const real_t, const len_t, MomentumGrid *mg, const RadialGrid*) {
    real_t
        *p    = new real_t[this->np],
        *p_f  = new real_t[this->np+1],
        *dp   = new real_t[this->np],
        *dp_f = nullptr;

    // Build flux grid
    if(pf_provided==nullptr)
        throw FVMException("PCustomGridGenerator: Rebuilding grid not yet supported.");
    for (len_t i = 0; i < this->np+1; i++)
        p_f[i] = this->pf_provided[i];
    pf_provided=nullptr;

    for (len_t i = 0; i < this->np; i++)
        dp[i] = p_f[i+1]-p_f[i];

    // Build cell grid
    for (len_t i = 0; i < this->np; i++)
        p[i] = 0.5 * (p_f[i+1] + p_f[i]);

    if (this->np > 1) {
        dp_f = new real_t[this->np-1];
        for (len_t i = 0; i < this->np-1; i++)
            dp_f[i] = p[i+1] - p[i];
    }

    mg->InitializeP1("p", this->np, p, p_f, dp, dp_f);

    initialized = true;
    return true;
}

