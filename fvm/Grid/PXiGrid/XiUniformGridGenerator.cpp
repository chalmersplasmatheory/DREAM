/**
 * Implementation of a simple uniform mometum (p) grid generator.
 */

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"


using namespace DREAM::FVM::PXiGrid;

/**
 * Constructor.
 *
 * nxi:   Number of points on (cell) grid.
 * xiMin: Minimum value of coordinate on (flux) grid.
 * xiMax: Maximum value of coordinate on (flux) grid.
 */
XiUniformGridGenerator::XiUniformGridGenerator(
    const len_t nxi
) : nxi(nxi) {
    
    if (nxi < 1)
        throw MomentumGridGeneratorException(
            "Uniform xi grid generator: A uniform grid must contain at least 1 cell. Specified number of cells: " LEN_T_PRINTF_FMT ".", nxi
        );

    if (xiMin == xiMax)
        throw MomentumGridGeneratorException(
            "Uniform xi grid generator: The two endpoints of a uniform grid may not be the equal."
        );

    if (xiMin > xiMax)
        throw MomentumGridGeneratorException(
            "Uniform xi grid generator: The lower boundary 'xiMin' of a uniform grid must be located below 'xiMax'."
        );
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
bool XiUniformGridGenerator::Rebuild(const real_t, const len_t, MomentumGrid *mg, const RadialGrid*) {
    real_t
        *xi    = new real_t[this->nxi],
        *xi_f  = new real_t[this->nxi+1],
        *dxi   = new real_t[this->nxi],
        *dxi_f = nullptr;

    // Build flux grid
    for (len_t i = 0; i < this->nxi; i++)
        dxi[i] = 2.0 / this->nxi;

    for (len_t i = 0; i < this->nxi+1; i++)
        xi_f[i] = (i * dxi[0]) - 1.0;

    // Build cell grid
    for (len_t i = 0; i < this->nxi; i++)
        xi[i] = 0.5 * (xi_f[i+1] + xi_f[i]);

    if (this->nxi > 1) {
        dxi_f = new real_t[this->nxi-1];
        for (len_t i = 0; i < this->nxi-1; i++)
            dxi_f[i] = xi[i+1] - xi[i];
    }

    mg->InitializeP2("xi", this->nxi, xi, xi_f, dxi, dxi_f);

    this->initialized = true;
    return true;
}

