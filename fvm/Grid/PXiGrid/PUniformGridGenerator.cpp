/**
 * Implementation of a simple uniform mometum (p) grid generator.
 */

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"


using namespace DREAM::FVM::PXiGrid;

/**
 * Constructor.
 *
 * np:   Number of points on (cell) grid.
 * pMin: Minimum value of coordinate on (flux) grid.
 * pMax: Maximum value of coordinate on (flux) grid.
 */
PUniformGridGenerator::PUniformGridGenerator(
    const len_t np, const real_t pMin, const real_t pMax
) : np(np), pMin(pMin), pMax(pMax) {
    
    if (np < 1)
        throw MomentumGridGeneratorException(
            "Uniform p grid generator: A uniform grid must contain at least 1 cell. Specified number of cells: " LEN_T_PRINTF_FMT ".", np
        );

    if (pMin == pMax)
        throw MomentumGridGeneratorException(
            "Uniform p grid generator: The two endpoints of a uniform grid may not be the equal."
        );

    if (pMin > pMax)
        throw MomentumGridGeneratorException(
            "Uniform p grid generator: The lower boundary 'pMin' of a uniform grid must be located below 'pMax'."
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
bool PUniformGridGenerator::Rebuild(const real_t, const len_t, MomentumGrid *mg, const RadialGrid*) {
    real_t
        *p    = new real_t[this->np],
        *p_f  = new real_t[this->np+1],
        *dp   = new real_t[this->np],
        *dp_f = nullptr;

    // Build flux grid
    for (len_t i = 0; i < this->np; i++)
        dp[i] = (this->pMax - this->pMin) / this->np;

    for (len_t i = 0; i < this->np+1; i++)
        p_f[i] = this->pMin + (i * dp[0]);

    // Build cell grid
    for (len_t i = 0; i < this->np; i++)
        p[i] = 0.5 * (p_f[i+1] + p_f[i]);

    if (this->np > 1) {
        dp_f = new real_t[this->np-1];
        for (len_t i = 0; i < this->np-1; i++)
            dp_f[i] = p[i+1] - p[i];
    }

    mg->InitializeP1(this->np, p, p_f, dp, dp_f);

    return true;
}

