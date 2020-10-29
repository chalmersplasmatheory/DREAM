/**
 * Implementation of a simple uniform mometum (p) grid generator.
 */

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiBiUniformThetaGridGenerator.hpp"
#include <cmath>

using namespace DREAM::FVM::PXiGrid;

/**
 * Constructor.
 *
 * nxi:         Number of points on (cell) grid.
 * nxiSep_frac: Fraction of grid points located below theta = thSep.
 * nxiSep:      Number of grid points located below theta = thSep.
 * thSep:       theta at which the grid spacing changes.
 */
XiBiUniformThetaGridGenerator::XiBiUniformThetaGridGenerator(
    const len_t nxi, const real_t nxiSep_frac, const real_t thSep
) : XiBiUniformThetaGridGenerator(nxi, static_cast<const len_t>(round(nxiSep_frac*nxi)), thSep) { }

XiBiUniformThetaGridGenerator::XiBiUniformThetaGridGenerator(
    const len_t nxi, const len_t nxiSep, const real_t thSep
) : nxi(nxi), nxiSep(nxiSep), thSep(thSep) {
    // Throw exceptions if invalid parameters.
    if ( (nxi < 2) || (nxiSep<1) )
        throw MomentumGridGeneratorException(
            "Two-region theta grid generator: A two-region grid must contain at least 1 cell in each region. Specified number of cells: " LEN_T_PRINTF_FMT " (lower) + " LEN_T_PRINTF_FMT " (total).", nxiSep, nxi
        );
    if (thMin == thMax)
        throw MomentumGridGeneratorException(
            "Two-region theta grid generator: The endpoints or interface of a two-region grid may not be the equal."
        );

    if (thMin > thMax)
        throw MomentumGridGeneratorException(
            "Two-region theta grid generator: The lower boundary 'xiMin' of a two-region grid must be located below 'xiMax'."
        );
    
    if (thSep > thMax)
        throw MomentumGridGeneratorException(
            "Two-region theta grid generator: The separating boundary 'thSep' of a two-region grid must be located below 'thMax'."
        );
    if (thMin > thSep)
        throw MomentumGridGeneratorException(
            "Two-region theta grid generator: The lower boundary 'thMin' of a two-region grid must be located below 'thSep'."
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
bool XiBiUniformThetaGridGenerator::Rebuild(const real_t, const len_t, MomentumGrid *mg, const RadialGrid*) {
    real_t
        *xi    = new real_t[this->nxi],
        *xi_f  = new real_t[this->nxi+1],
        *th_f  = new real_t[this->nxi+1],
        *dxi   = new real_t[this->nxi],
        *dth   = new real_t[this->nxi],
        *dxi_f = nullptr;
        
    real_t dth_lo = (thSep - thMin) / this->nxiSep;
    real_t dth_up  = (thMax - thSep) / (this->nxi - this->nxiSep);

    // Build flux grid
    for (len_t i = 0; i < this->nxiSep; i++)
        dth[i] = dth_lo;
    for (len_t i = this->nxiSep; i < this->nxi; i++)
        dth[i] = dth_up;
	
    th_f[0] = thMin;
    xi_f[this->nxi] = cos(thMin);
    xi_f[0] = cos(thMax);
   
    for (len_t i = 1; i <= this->nxi-1; i++){
        th_f[i] = th_f[i-1] + dth[i-1];
        xi_f[this->nxi-i] = cos(th_f[i]);
        }  
    
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
    delete[]dth;
    delete[]th_f;
    this->initialized = true;
    return true;
}

