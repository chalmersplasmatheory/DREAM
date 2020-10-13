/**
 * Implementation of the Connor-Hastie runaway rate.
 */

#include "FVM/config.h"
#include "DREAM/Equations/ConnorHastie.hpp"

using namespace DREAM;


/**
 * Constructor.
 *
 * rf:          Parent RunawayFluid object.
 * corrections: If 'true', includes corrections for the parameters
 *              'lambda', 'eta' and 'h' appearing in the runaway rate.
 */
ConnorHastie::ConnorHastie(RunawayFluid *rf, bool corrections)
    : REFluid(rf), withCorrections(corrections) { }


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
    if (Ec >= abs(E)) return 0;

    real_t ED = REFluid->GetDreicerElectricField(ir);
    real_t tauEE = REFluid->GetElectronCollisionTimeThermal(ir);
    
    real_t EEc   = abs(E) / Ec;
    real_t EED   = abs(E) / ED;
    
    // "Undetermined" factor (~1 is usually good according to simulations)
    real_t C = 1;
    real_t h=1, eta=1, etaf=1, lmbd=1;

    // Include corrections?
    if (this->withCorrections) {
        h    = 1.0/(3*(EEc-1)) * (EEc + 2*(EEc-2)*sqrt(EEc/(EEc-1)) - (Zeff-7)/(Zeff+1));
        etaf = (M_PI/2 - asin(1-2/EEc));
        eta  = EEc*EEc/(4*(EEc-1)) * etaf*etaf;
        lmbd = 8*EEc*EEc*(1 - 1.0/(2*EEc) - sqrt(1-1/EEc));
    }
    
    real_t alpha = -3.0/16.0*(1+Zeff)*h;
    real_t baseExpr =
        C*ne/tauEE * pow(EED, alpha) *
        exp(-lmbd/(4*EED) - sqrt(eta*(1+Zeff)/EED));

    if (derivative)
        return baseExpr *
            (alpha/EED + lmbd/(4*EED*EED) + 0.5*sqrt(eta*(1+Zeff)/EED)/EED);
    else
        return baseExpr;
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

