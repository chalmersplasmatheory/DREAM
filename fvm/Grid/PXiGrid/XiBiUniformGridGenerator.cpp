/**
 * Implementation of a simple uniform mometum (p) grid generator.
 */

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiBiUniformGridGenerator.hpp"


using namespace DREAM::FVM::PXiGrid;

/**
 * Constructor.
 *
 * nxi:   Number of points on (cell) grid.
 * xiMin: Minimum value of coordinate on (flux) grid.
 * xiMax: Maximum value of coordinate on (flux) grid.
 */
XiBiUniformGridGenerator::XiBiUniformGridGenerator(
    const len_t nxi, const len_t nxiSep, const real_t xiSep
) : nxi(nxi), nxiSep(nxiSep), xiSep(xiSep) {
    // Throw exceptions if invalid parameters.
    if ( (nxi < 2) || (nxiSep<1) )
        throw MomentumGridGeneratorException(
            "Two-region xi grid generator: A two-region grid must contain at least 1 cell in each region. Specified number of cells: " LEN_T_PRINTF_FMT " (lower) + " LEN_T_PRINTF_FMT " (total).", nxiSep, nxi
        );
    if (xiMin == xiMax)
        throw MomentumGridGeneratorException(
            "Two-region xi grid generator: The endpoints or interface of a two-region grid may not be the equal."
        );

    if (xiMin > xiMax)
        throw MomentumGridGeneratorException(
            "Two-region xi grid generator: The lower boundary 'xiMin' of a two-region grid must be located below 'xiMax'."
        );
    
    if (xiSep > xiMax)
        throw MomentumGridGeneratorException(
            "Two-region xi grid generator: The separating boundary 'xiSep' of a two-region grid must be located below 'xiMax'."
        );
    if (xiMin > xiSep)
        throw MomentumGridGeneratorException(
            "Two-region xi grid generator: The lower boundary 'xiMin' of a two-region grid must be located below 'xiSep'."
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
bool XiBiUniformGridGenerator::Rebuild(const real_t, const len_t, MomentumGrid *mg, const RadialGrid*) {
    real_t
        *xi    = new real_t[this->nxi],
        *xi_f  = new real_t[this->nxi+1],
        *dxi   = new real_t[this->nxi],
        *dxi_f = nullptr;

    real_t dxi_lo = (this->xiSep - this->xiMin) / this->nxiSep;
    real_t dxi_up  = (this->xiMax - this->xiSep) / (this->nxi - this->nxiSep);


    // Build flux grid
    for (len_t i = 0; i < this->nxiSep; i++)
        dxi[i] = dxi_lo;
    for (len_t i = this->nxiSep; i < this->nxi; i++)
        dxi[i] = dxi_up;

    xi_f[0] = this->xiMin;
    for (len_t i = 1; i <= this->nxi; i++)
        xi_f[i] = xi_f[i-1] + dxi[i-1];

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

