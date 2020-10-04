/**
 * This equation term implements turbulent electron radial diffusion as
 * described by Rechester & Rosenbluth (PRL 40 1978):
 *
 *   Drr = pi*q*R0*(dB/B)^2 * v_||
 *
 * where q = q(r) denotes the tokamak safety factor, R0 is the tokamak major
 * radius, v_|| is the electron parallel speed, and dB is the magnitude by which
 * the magnetic field strength oscillates.
 */


#include "DREAM/Equations/Kinetic/RechesterRosenbluthTransport.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;


/**
 * Constructor.
 *
 * grid: Computational grid on which the equation term is defined.
 * dB_B: Object representing the time evolution of the magntiude of
 *       the magnetic field oscillation (dB/B).
 */
RechesterRosenbluthTransport::RechesterRosenbluthTransport(
    FVM::Grid *grid, enum OptionConstants::momentumgrid_type mgtype, FVM::Interpolator1D *dB_B
) : DiffusionTerm(grid), mgtype(mgtype), deltaBOverB(dB_B) { }

/**
 * Destructor.
 */
RechesterRosenbluthTransport::~RechesterRosenbluthTransport() {
    delete this->deltaBOverB;
}


/**
 * Rebuild this equation term.
 */
void RechesterRosenbluthTransport::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler*
) {
    const real_t *dB_B = this->deltaBOverB->Eval(t);

    // XXX Here we assume that all momentum grids are the same
    // at all radii...
    const len_t nr = this->grid->GetNr();
    auto mg = this->grid->GetMomentumGrid(0);

    const real_t
        *p1 = mg->GetP1(),
        *p2 = mg->GetP2(),
        R0  = this->grid->GetRadialGrid()->GetR0();
    const len_t
        np1 = mg->GetNp1(),
        np2 = mg->GetNp2();

    for (len_t ir = 0; ir < nr; ir++) {
        const real_t q = 1.0;   // TODO (safety factor)

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                real_t vpar=1;
                if (mgtype == OptionConstants::MOMENTUMGRID_TYPE_PXI)
                    vpar = p1[i]*p2[j] / sqrt(1+p1[i]*p1[i]);
                else if (mgtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP)
                    vpar = p1[i] / sqrt(1 + p1[i]*p1[i] + p2[j]*p2[j]);

                // Set diffusion coefficient...
                Drr(ir, i, j) = M_PI * q * R0 * dB_B[ir]*dB_B[ir] * vpar;
            }
        }
    }
}

