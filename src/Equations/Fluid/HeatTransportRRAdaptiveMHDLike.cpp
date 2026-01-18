/**
 * Adaptive MHD-like heat transport, using a Rechester-Rosenbluth
 * diffusion coefficient.
 */

#include "DREAM/Equations/Fluid/HeatTransportRRAdaptiveMHDLike.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
HeatTransportRRAdaptiveMHDLike::HeatTransportRRAdaptiveMHDLike(
    FVM::Grid *grid,
    FVM::UnknownQuantityHandler *uqh,
    RunawayFluid *REFluid,
    real_t L0, real_t Lperp,
    const real_t grad_j_tot_max,
    bool gradient_normalized,
    const real_t dBOverB,
    const real_t suppression_level,
    bool localized
) : AdaptiveMHDLikeTransportTerm(
        grid, uqh,
        grad_j_tot_max, gradient_normalized,
        suppression_level, localized
    ),
    // === NEW ===============================================================
    // Passing REFluid enables the RR base to access the thermal collision time
    // tau_th(r) and build nu_c(r)=1/tau_th(r) consistently.
    // ======================================================================
    HeatTransportRechesterRosenbluth(
        grid,
        nullptr,      // dB_B interpolator (unused, overridden by EvaluateDeltaBOverB)
        uqh,
        REFluid,
        L0, Lperp
    ),
    dBOverB(dBOverB)
{
    this->dB = new real_t[grid->GetNr()];
}


/**
 * Destructor.
 */
HeatTransportRRAdaptiveMHDLike::~HeatTransportRRAdaptiveMHDLike() {
    if (this->dB != nullptr)
        delete [] this->dB;
}


/**
 * Evaluate dB/B.
 */
const real_t *HeatTransportRRAdaptiveMHDLike::EvaluateDeltaBOverB(const real_t t) {

    // === NEW ===============================================================
    // In the adaptive model, transport is enabled only when the trigger
    // condition is satisfied (e.g. when |∇j| exceeds a threshold). If transport
    // is disabled, δB/B is set to zero everywhere, effectively turning off RR
    // diffusion in the base class.
    // ======================================================================

    real_t v = 0;
    if (this->CheckTransportEnabled(t))
        v = this->dBOverB;

    // === NEW ===============================================================
    // Apply a spatial mask to localize transport. The mask[] array is managed
    // by AdaptiveMHDLikeTransportTerm ( selecting regions where enhanced
    // transport should be active). The returned δB/B profile is:
    //     (δB/B)(t,r) = v(t) * mask(r)
    // and is then used by HeatTransportRechesterRosenbluth::Rebuild().
    // ======================================================================
    
    const len_t nr = this->AdaptiveMHDLikeTransportTerm::grid->GetNr();
    for (len_t ir = 0; ir < nr; ir++)
        this->dB[ir] = v * this->mask[ir];
    
    return this->dB;
}
