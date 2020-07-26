/**
 * Implementation of the Connor-Hastie runaway rate.
 */

#include "FVM/config.h"
#include "DREAM/Equations/ConnorHastie.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
ConnorHastie::ConnorHastie(RunawayFluid *rf)
    : REFluid(rf) { }


/**
 * Evaluate the Connor-Hastie runaway rate.
 *
 * ir:    Index of radius for which to evaluate the runaway rate.
 * E:     Local electric field strength.
 * ne:    Electron density.
 * Zeff:  Local effective plasma charge.
 * deriv: If true, returns the runaway rate derivative with respect
 *        to E/ED (electric field divided by Dreicer field) instead
 *        of the runaway rate.
 */
real_t ConnorHastie::RunawayRate(
    const len_t ir, const real_t E, const real_t ne, const real_t Zeff,
    bool derivative
) {
    if (E == 0) return 0; 
    
    real_t Ec = REFluid->GetConnorHastieField_COMPLETESCREENING(ir);
    if (Ec >= E) return 0;

    real_t ED = REFluid->GetDreicerElectricField(ir);
    real_t tauEE = REFluid->GetElectronCollisionTimeThermal(ir);
    
    real_t EEc   = E / Ec;
    real_t EED   = E / ED;
    
    // "Undetermined" factor (~1 is usually good according to simulations)
    real_t C = 1;

    real_t h    = 1.0/(3*(EEc-1)) * (EEc + 2*(EEc-2)*sqrt(EEc/(EEc-1)) - (Zeff-7)/(Zeff+1));
    real_t etaf = (M_PI/2 - asin(1-2/EEc));
    real_t eta  = EEc*EEc/(4*(EEc-1)) * etaf*etaf;
    real_t lmbd = 8*EEc*EEc*(1 - 1.0/(2*EEc) - sqrt(1-1/EEc));
    
    real_t alpha = 3.0/16.0*(1+Zeff)*h;

    if (derivative)
        return C*ne/tauEE * pow(EED, alpha) *
            (alpha/EED - lmbd/(4*EED*EED) - 0.5*sqrt(eta*(1+Zeff)/EED)/EED) *
            exp(-lmbd/(4*EED) - sqrt(eta*(1+Zeff)/EED));
    else
        return C*ne/tauEE * pow(EED, alpha) *
            exp(-lmbd/(4*EED) - sqrt(eta*(1+Zeff)/EED));
}

/**
 * Derivative of runaway rate with respect to E/E_D.
 *
 *   d (gamma) / d (E/E_D)
 *
 * ir:   Index of radius for which to evaluate the runaway rate.
 * E:    Local electric field strength.
 * ne:   Electron density.
 * Zeff: Local effective plasma charge.
 */
real_t ConnorHastie::Diff_EED(
    const len_t ir, const real_t E, const real_t ne, const real_t Zeff
) {
    return RunawayRate(ir, E, ne, Zeff, true);
}

/**
 * Derivative of runaway rate with respect to E.
 *
 * ir:   Index of radius for which to evaluate the runaway rate.
 * E:    Local electric field strength.
 * ne:   Electron density.
 * Zeff: Local effective plasma charge.
 */
real_t ConnorHastie::Diff_E(
    const len_t ir, const real_t E, const real_t ne, const real_t Zeff
) {
    real_t ED = REFluid->GetDreicerElectricField(ir);

    return Diff_EED(ir, E, ne, Zeff) / ED;
}

/**
 * Derivative of runaway rate with respect to ne.
 *
 * ir:   Index of radius for which to evaluate the runaway rate.
 * E:    Local electric field strength.
 * ne:   Electron density.
 * Zeff: Local effective plasma charge.
 */
real_t ConnorHastie::Diff_ne(
    const len_t ir, const real_t E, const real_t ne, const real_t Zeff
) {
    real_t ED = REFluid->GetDreicerElectricField(ir);

    return Diff_EED(ir, E, ne, Zeff) * (E/ED) / ne;
}

/**
 * Derivative of runaway rate with respect to temperature Te.
 *
 * ir:   Index of radius for which to evaluate the runaway rate.
 * E:    Local electric field strength.
 * ne:   Electron density.
 * Zeff: Local effective plasma charge.
 * Te:   Electron temperature.
 */
real_t ConnorHastie::Diff_Te(
    const len_t ir, const real_t E, const real_t ne, const real_t Zeff,
    const real_t Te
) {
    real_t ED = REFluid->GetDreicerElectricField(ir);

    return Diff_EED(ir, E, ne, Zeff) * (E/ED) / Te;
}

