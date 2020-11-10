/**
 * Implementaion of base class from which fluid runaway source terms
 * should be derived. This class contains routines which can be useful
 * to most or all of the various runaway source terms.
 */

#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
RunawaySourceTerm::RunawaySourceTerm(
    FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns
) : rst_grid(grid), rst_unknowns(unknowns) {
    this->id_E_field = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
}

/**
 * Returns the xi index for either xi=-1 or xi=+1, depending on the sign of
 * the electric field at the given radius.
 *
 * ir: Radial index.
 */
len_t RunawaySourceTerm::GetXiIndexForEDirection(const len_t ir) {
    real_t E = this->rst_unknowns->GetUnknownData(this->id_E_field)[ir];

    if (E < 0)
        return 0;
    else
        return this->rst_grid->GetMomentumGrid(ir)->GetNp2()-1;
}

/**
 * Returns the appropriate volume element scale factor to use for the
 * term. If this term is applied to a fluid grid, then the scale factor
 * is 1 (no re-scaling), but if it's applied to a kinetic grid, the scale
 * factor is "VpVol/(dxi*dp*Vp)" in order for the correct number of
 * particles to be added to the grid.
 *
 * On a kinetic grid, the runaway source should be placed fully inside the
 * cell at p_{1/2} and xi=+1 (if E>0) or xi=-1 (if E<0).
 */
real_t RunawaySourceTerm::GetVolumeScaleFactor(const len_t ir) {
    const len_t np1 = this->rst_grid->GetMomentumGrid(ir)->GetNp1();
    const len_t np2 = this->rst_grid->GetMomentumGrid(ir)->GetNp2();

    const len_t xiIndex = this->GetXiIndexForEDirection(ir);

    // Kinetic grid?
    if (np1 == 1 && np2 == 1)
        return 1.0;

    real_t VpVol = this->rst_grid->GetVpVol(ir);
    real_t Vp    = this->rst_grid->GetVp(ir, 0, xiIndex);
    real_t dp    = this->rst_grid->GetMomentumGrid(ir)->GetDp1(0);
    real_t dxi   = this->rst_grid->GetMomentumGrid(ir)->GetDp2(xiIndex);

    return VpVol / (dxi*dp*Vp);
}

