/**
 * Implementation of a bi-uniform momentum (p) grid generator.
 * It describes a piecewise constant with two regions
 * separated by momentum pSep.
 */

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PBiUniformGridGenerator.hpp"


using namespace DREAM::FVM::PXiGrid;

/**
 * Constructor.
 *
 * np:         Number of points on (cell) grid.
 * npSep_frac: Fraction of grid points below p = pSep.
 * npSep:      Number of grid points below p = pSep.
 * pMin:       Minimum value of coordinate on (flux) grid.
 * pSep:       Grid point at which the grid spacing changes.
 * pMax:       Maximum value of coordinate on (flux) grid.
 */
PBiUniformGridGenerator::PBiUniformGridGenerator(
    const len_t np, const real_t npSep_frac, const real_t pMin, const real_t pSep, const real_t pMax
) : PBiUniformGridGenerator(np, static_cast<const len_t>(round(npSep_frac*np)), pMin, pSep, pMax) { }

PBiUniformGridGenerator::PBiUniformGridGenerator(
    const len_t np, const len_t npSep, const real_t pMin, const real_t pSep, const real_t pMax
) :  np(np), npSep(npSep), pMin(pMin), pSep(pSep), pMax(pMax) {
    // Throw exceptions if invalid parameters.
    if ( (np < 2) || (npSep<1) )
        throw MomentumGridGeneratorException(
            "Two-region p grid generator: A two-region grid must contain at least 1 cell in each region. Specified number of cells: " LEN_T_PRINTF_FMT " (lower) + " LEN_T_PRINTF_FMT " (total).", npSep, np
        );
    if (( pMin == pMax) || (pMin == pSep) || (pSep == pMax) ) 
        throw MomentumGridGeneratorException(
            "Two-region p grid generator: The endpoints or interface of a two-region grid may not be the equal."
        );
    if (pMin > pMax)
        throw MomentumGridGeneratorException(
            "Two-region p grid generator: The lower boundary 'pMin' of a two-region grid must be located below 'pMax'."
        );
    if (pSep > pMax)
        throw MomentumGridGeneratorException(
            "Two-region p grid generator: The separating boundary 'pSep' of a two-region grid must be located below 'pMax'."
        );
    if (pMin > pSep)
        throw MomentumGridGeneratorException(
            "Two-region p grid generator: The lower boundary 'pMin' of a two-region grid must be located below 'pSep'."
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
bool PBiUniformGridGenerator::Rebuild(const real_t, const len_t, MomentumGrid *mg, const RadialGrid*) {
    real_t
        *p    = new real_t[this->np],
        *p_f  = new real_t[this->np+1],
        *dp   = new real_t[this->np],
        *dp_f = nullptr;

    real_t dp_lo = (this->pSep - this->pMin) / this->npSep;
    real_t dp_up  = (this->pMax - this->pSep) / (this->np - this->npSep);
    
    // Build flux grid
    for (len_t i = 0; i < this->npSep; i++)
        dp[i] = dp_lo;
    for (len_t i = this->npSep; i < this->np; i++)
        dp[i] = dp_up;

    p_f[0] = this->pMin;
    for (len_t i = 1; i < this->np+1; i++)
        p_f[i] = p_f[i-1] + dp[i-1];

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

