/**
 * Implementation of a class that calculates and stores quantities related to
 * bootstrap current. For further details on the implementation, see doc/notes/bootstrap.
 */
#include "DREAM/Equations/BootstrapCurrent.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Constants.hpp"
#include <iostream>
using namespace DREAM;

/**
 * Constructor.
 */
BootstrapCurrent::BootstrapCurrent(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih, CoulombLogarithm *lnL, enum OptionConstants::eqterm_bootstrap_mode mode) {

    rGrid = g->GetRadialGrid();
    unknowns = u;
    ions = ih;
    lnLambda = lnL;


    // used for partial derivatives of lnLambda
    lnLII_settings = new CollisionQuantity::collqty_settings;
    lnLII_settings->lnL_type = OptionConstants::COLLQTY_LNLAMBDA_ION_ION;
    lnLEE_settings = new CollisionQuantity::collqty_settings;
    lnLEE_settings->lnL_type = OptionConstants::COLLQTY_LNLAMBDA_THERMAL;

    // cube root of a small number used in central difference derivatives
    epsilon_forward = sqrt(std::numeric_limits<real_t>::epsilon());
    epsilon_central = cbrt(std::numeric_limits<real_t>::epsilon());

    id_jtot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ions  = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);   // cc summed densities
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

    includeIonTemperatures = u->HasUnknown(OptionConstants::UQTY_WI_ENER);
    if (includeIonTemperatures)
        id_Wi = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);

    nr = rGrid->GetNr();

    // Memory allocation
    AllocateQuantities();

    // equilibrium constants
    const real_t R0 = rGrid->GetR0();
    const real_t B0 = rGrid->GetBmin_f(0);
    if (mode == OptionConstants::EQTERM_BOOTSTRAP_MODE_REDL) {
        for (len_t ir = 0; ir < nr; ir++) {
            // calculate the geometric prefactor
            const real_t BtorGOverR0 = rGrid->GetBTorG(ir);        // G / R0
            const real_t FSA_B2 = rGrid->GetFSA_B2(ir);            // <B^2> / Bmin^2
            const real_t Bmin = rGrid->GetBmin(ir);                // Bmin
            const real_t psiPrimeRef = rGrid->GetPsiPrimeRef(ir);  // R0 d(psi_ref)/dr

            // OBS. something is off with the above definitions: the following should not include the last factor of B0...
            constantPrefactor[ir] = -BtorGOverR0 * R0 * R0 / ( FSA_B2 * Bmin * psiPrimeRef)  *B0; // <--- this B0!
            if (ir == 0)
                constantPrefactor[ir] /= 2 * rGrid->GetDr_f(ir);
            else if (ir == nr - 1)
                constantPrefactor[ir] /= rGrid->GetDr_f(ir-1);
            else
                constantPrefactor[ir] /= ( rGrid->GetDr_f(ir-1) + rGrid->GetDr_f(ir) );

            // convert eV to J by multiplying with the electron charge
            constantPrefactor[ir] *= Constants::ec;

            // calculate fraction of trapped particles
            ft[ir] = 1. - rGrid->GetEffPassFrac(ir);
            // ft[ir] = 1.46 * sqrt( rGrid->GetR(ir) / R0);

            // this high-aspect ratio approximation for qR0 seems to match better with Redl-Sauter
            // than calculating it via the total current (also simpler for initialization)
            qR0[ir] = rGrid->GetR(ir) * sqrt(1 + 4*M_PI*M_PI * BtorGOverR0 * BtorGOverR0 / (psiPrimeRef * psiPrimeRef));

            eps[ir] = rGrid->GetR(ir) / rGrid->GetR0();
        }
    } else if (mode == OptionConstants::EQTERM_BOOTSTRAP_MODE_REDL_STELLARATOR) {
        for (len_t ir = 0; ir < nr; ir++) {
            // calculate the geometric prefactor
            const real_t BtorGOverR0 = rGrid->GetBTorG(ir);        // G / R0
            const real_t FSA_B2 = rGrid->GetFSA_B2(ir);            // <B^2> / Bmin^2
            const real_t FSA_1OverB = rGrid->GetFSA_1OverB(ir);    // <1 / B> * Bmin
            const real_t Bmin = rGrid->GetBmin(ir);                // Bmin
            const real_t psiPrimeRef = rGrid->GetPsiPrimeRef(ir);  // R0 d(psi_ref)/dr

            // OBS. something is off with the above definitions: the following should not include the last factor of B0...
            constantPrefactor[ir] = -BtorGOverR0 * R0 * R0 / ( FSA_B2 * Bmin * psiPrimeRef)  *B0; // <--- this B0!
            if (ir == 0)
                constantPrefactor[ir] /= 2 * rGrid->GetDr_f(ir);
            else if (ir == nr - 1)
                constantPrefactor[ir] /= rGrid->GetDr_f(ir-1);
            else
                constantPrefactor[ir] /= ( rGrid->GetDr_f(ir-1) + rGrid->GetDr_f(ir) );

            // convert eV to J by multiplying with the electron charge
            constantPrefactor[ir] *= Constants::ec;

            // For stellarators, density and temperature gradients dX/dr->(dX/dr)/iota in the Redl formula
            constantPrefactor[ir] *= rGrid->GetIota(ir);

            // calculate fraction of trapped particles
            ft[ir] = 1. - rGrid->GetEffPassFrac(ir);
            // ft[ir] = 1.46 * sqrt( rGrid->GetR(ir) / R0);

            qR0[ir] = BtorGOverR0 * R0 / rGrid->GetIota(ir) * FSA_1OverB / Bmin; // IE: Should we divide by Bmin here?

            eps[ir] = (rGrid->GetBmax(ir) - rGrid->GetBmin(ir)) / (rGrid->GetBmax(ir) + rGrid->GetBmin(ir));
        }
    } /*else {
        Possibly do a warning here?
    }*/
    // IE: Do we want something else for a stellarator?

    // locate the main ion index
    bool isFound = false;
    // IE: is this (below) the best practice?
    for (len_t iZ = 0; iZ < ions->GetNZ(); iZ++)
        if (ions->GetZ(iZ) == 1) {
            iZMain = iZ;
            isFound = true;
            break;
        }
    if (!isFound)
        throw DREAMException("No hydrogenic ion was found to set as the main ion for calculating the bootstrap current!");
}

/**
 * Destructor.
 */
BootstrapCurrent::~BootstrapCurrent() {
    DeallocateQuantities();
    delete lnLII_settings;
    delete lnLEE_settings;
}

/**
 * Allocate memory for arrays stored in this object.
 */
void BootstrapCurrent::AllocateQuantities() {
    constantPrefactor = new real_t[nr];
    coefficientL31    = new real_t[nr];
    coefficientL32    = new real_t[nr];
    coefficientAlpha  = new real_t[nr];
    NiMain            = new real_t[nr];
    WiMain            = new real_t[nr];
    qR0               = new real_t[nr];
    eps               = new real_t[nr];
    ft                = new real_t[nr];
    n                 = new real_t[nr];
    p                 = new real_t[nr];
}

/**
 * Deallocate memory for arrays stored in this object.
 */
void BootstrapCurrent::DeallocateQuantities() {
    delete [] constantPrefactor;
    delete [] coefficientL31;
    delete [] coefficientL32;
    delete [] coefficientAlpha;
    delete [] NiMain;
    delete [] WiMain;
    delete [] qR0;
    delete [] eps;
    delete [] ft;
    delete [] n;
    delete [] p;
}

/**
 * Rebuild bootstrap coefficients and other relevant quantities.
 */
void BootstrapCurrent::Rebuild() {

    // jtot   = unknowns->GetUnknownData(id_jtot); // dependency neglected in jacobian
    ncold  = unknowns->GetUnknownData(id_ncold);
    Ni     = unknowns->GetUnknownData(id_Ni);
    Tcold  = unknowns->GetUnknownData(id_Tcold);
    if (includeIonTemperatures)
        Wi = unknowns->GetUnknownData(id_Wi);

    // check for when jtot is initialised, then use it to obtain the safety factor
    // if (!qFromCurrent)
    //     for (len_t ir = 0; ir < nr; ir++)
    //         if (jtot[ir] != 0)
    //             qFromCurrent = true;

    for (len_t ir = 0; ir < nr; ir++) {

        NiMain[ir] = Ni[nr * iZMain + ir];  // main ion density (for alpha)
        if (includeIonTemperatures)
            WiMain[ir] = Wi[nr * iZMain + ir] / Constants::ec; // main ion thermal energy (for alpha)
        else
            WiMain[ir] = 1.5 * NiMain[ir] * Tcold[ir];

        // Calculate number and pressure density
        n[ir] = ncold[ir];
        p[ir] = ncold[ir] * Tcold[ir];
        for (len_t i = ir; i < nr * ions->GetNZ(); i += nr) {
            n[ir] += Ni[i];
            if (includeIonTemperatures)
                p[ir] += Wi[i] / (1.5 * Constants::ec);
            else
                p[ir] += Ni[i] * Tcold[ir];
        }

        // if (qFromCurrent) {
        //     real_t mu0Ip = Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir, rGrid, jtot);
        //     qR0[ir] = fabs(rGrid->SafetyFactorNormalized(ir, mu0Ip));
        // }

        // calculate collision frequencies
        real_t nuE = evaluateElectronCollisionFrequency(ir);
        real_t nuI = evaluateIonCollisionFrequency(ir);

        // get the effective charge
        real_t Zeff =  ions->GetZeff(ir);

        // calculate the bootstrap coefficients
        coefficientL31[ir] = evaluateCoefficientL31(ft[ir], Zeff, nuE);
        coefficientL32[ir] = evaluateCoefficientL32(ft[ir], Zeff, nuE);
        coefficientAlpha[ir] = evaluateCoefficientAlpha(ft[ir], Zeff, nuI);
    }
}


/**
 * Thermal electron collision frequency (as defined in Eq. 18b(d) in Sauter et al. 1999).
 *
 * ir:  radial cell grid point.
 */
real_t BootstrapCurrent::evaluateElectronCollisionFrequency(len_t ir) {
    real_t lnLee = lnLambda->evaluateLnLambdaT(ir);
    real_t Zeff = ions->GetZeff(ir);
    return 6.921e-18 * ncold[ir] * lnLee * Zeff * qR0[ir] / (eps[ir] * sqrt(eps[ir]) * Tcold[ir] * Tcold[ir]);
}

/**
 * Ion collision frequency ( as defined in Eq. 18c(e) in Sauter et al. 1999).
 *
 * ir:  radial cell grid point.
 */
real_t BootstrapCurrent::evaluateIonCollisionFrequency(len_t ir) {
    real_t lnLii = lnLambda->evaluateLnLambdaII(ir); // formula from Wesson
    real_t TiMain = WiMain[ir] / (1.5 * NiMain[ir]);
    real_t TiMain2 = TiMain * TiMain;
    real_t Zeff = ions->GetZeff(ir);
    real_t Zeff4 = Zeff * Zeff * Zeff * Zeff;
    return 4.90e-18  * NiMain[ir] * lnLii * Zeff4 * qR0[ir] / (eps[ir] * sqrt(eps[ir]) * TiMain2 );
}


//
// METHODS FOR CALCULATING THE BOOTSTRAP COEFFICIENTS
//

/**
 * Calculates the coefficient L31 (as defined in Eqs. 10-11 in Redl et al. 2021).
 *
 * ft:    fraction of trapped particles.
 * Zeff:  effective ion charge.
 * nu:    electron collision frequency.
 */
real_t BootstrapCurrent::evaluateCoefficientL31(real_t ft, real_t Zeff, real_t nu) {

    // Eq. 11
    real_t d = (.52 + .086 * sqrt(nu)) * (1. + .87 * ft) * nu / (1. + 1.13 * sqrt(Zeff - 1.));
    d += 1 + .67*(1. - .7 * ft) * sqrt(nu) / (.56 + .44*Zeff);
    real_t f31 = ft / d;

    // Eq. 10
    real_t a = 1 / (pow(Zeff, 1.2) - .71);
    return (1. + .15*a)*f31 - .22*a*f31*f31 + .01*a*f31*f31*f31 + .06*a*f31*f31*f31*f31;
}


/**
 * Calculates the coefficient L32 (as defined in Eqs. 12-16 in Redl et al. 2021).
 *
 * ft:    fraction of trapped particles.
 * Zeff:  effective ion charge.
 * nu:    electron collision frequency.
 */
real_t BootstrapCurrent::evaluateCoefficientL32(real_t ft, real_t Zeff, real_t nu) {

    // Eq. 14
    real_t dee = sqrt( 1. + 2. * sqrt(Zeff - 1.) );
    dee += ft * ft * sqrt( nu * (.075 + .25 * (Zeff - 1.) * ( Zeff - 1.)) );
    dee *= .13 * (1. - .38 * ft) * nu / (Zeff * Zeff);
    dee += 1. + .23 * (1. - .96 * ft) * sqrt(nu / Zeff);
    real_t f32ee = ft / dee;

    // Eq. 13
    real_t F32ee = (.1 + .6 * Zeff) * f32ee * (1. - f32ee * f32ee * f32ee);
    F32ee /= Zeff * (.77 + .63*(1. + pow(Zeff-1., 1.1)));
    F32ee += .7 / (1. + .2 * Zeff) * f32ee * f32ee * (1. - 1.2 * f32ee + .2 * f32ee * f32ee);
    F32ee += 1.3 * f32ee * f32ee * f32ee * f32ee / (1. + .5 * Zeff);

    // Eq. 16
    real_t dei = 1. + .87 * (1. + .39 * ft) * sqrt(nu) / (1. + 2.95 * (Zeff - 1.) * (Zeff - 1.));
    dei += 1.53 * (1. - .37 * ft) * nu * (2. + .375 * (Zeff - 1.));
    real_t f32ei = ft / dei;

    // Eq. 15
    real_t F32ei = 5.5 * f32ei * f32ei * (1. - .8 * f32ei - .2 * f32ei * f32ei);
    F32ei /= (1.5 + 2. * Zeff);
    F32ei -= (.4 + 1.93 * Zeff) * f32ei * (1. - f32ei * f32ei * f32ei) / (Zeff * (.8 + .6 * Zeff));
    F32ei -= 1.3 * f32ei * f32ei * f32ei * f32ei / (1. + .5 * Zeff);

    return F32ee + F32ei;
}

/**
 * Calculates the coefficient alpha (as defined in Eqs. 20-21 in Redl et al. 2021).
 *
 * ft:    fraction of trapped particles.
 * Zeff:  effective ion charge.
 * nu:    ion collision frequency.
 */
real_t BootstrapCurrent::evaluateCoefficientAlpha(real_t ft, real_t Zeff, real_t nu) {

    // Eq. 20
    real_t alpha0 = - ( .62 + .055 * (Zeff - 1.) ) / ( .53 + .17 * (Zeff - 1.) );
    real_t ft2 = ft * ft;
    alpha0 *= (1. - ft) / (1. - ( .31 - .065 * (Zeff - 1.)) * ft - .25 * ft2);

    // Eq. 21
    real_t alpha = ( alpha0 + .7 * Zeff * sqrt(ft * nu) ) / ( 1. + .18 * sqrt(nu) );
    real_t ft6 = ft2 * ft2 * ft2;
    alpha -= .002 * nu*nu * ft6;
    alpha /= 1. + .004 * nu*nu * ft6;
    return alpha;
}




//
// METHODS FOR CALCULATING PARTIAL DERIVATIVES OF THE BOOTSTRAP COEFFICIENTS
//


/**
 * Calculates the partial derivate of the coefficient L31.
 *
 *  ir:       radial cell grid point.
 *  derivId:  unknown quantity ID to differentiate with respect to.
 *  iz:       ion charge state index.
 */
real_t BootstrapCurrent::evaluatePartialCoefficientL31(len_t ir, len_t derivId, len_t iz) {
    return evaluateNumericalDerivative(ir, derivId, iz, &evaluateCoefficientL31);
}


/**
 * Calculates the partial derivate of the coefficient L32.
 *
 *  ir:       radial cell grid point.
 *  derivId:  unknown quantity ID to differentiate with respect to.
 *  iz:       ion charge state index.
 */
real_t BootstrapCurrent::evaluatePartialCoefficientL32(len_t ir, len_t derivId, len_t iz) {
    return evaluateNumericalDerivative(ir, derivId, iz, &evaluateCoefficientL32);
}


/**
 * Numerical differentiation used for evaluating partial derivatives of coefficient L31 or L32.
 *
 *  ir:           radial cell grid point.
 *  derivId:      unknown quantity ID to differentiate with respect to.
 *  iz:           ion charge state index.
 *  coefficient:  coefficient function of variables (ir, Zeff, nu) to differentiate.
 */
real_t BootstrapCurrent::evaluateNumericalDerivative(
    len_t ir, len_t derivId, len_t izs,
    std::function<real_t(real_t, real_t, real_t)> coefficient
) {
    if ( (derivId != id_ncold) && (derivId != id_Tcold) && (derivId != id_ions) )
        return 0;
    real_t nu = evaluateElectronCollisionFrequency(ir);
    real_t Zeff = ions->GetZeff(ir);
    real_t hnu = nu * epsilon_central;
    real_t dCdnu = ( coefficient(ft[ir], Zeff, nu+hnu)
                   - coefficient(ft[ir], Zeff, nu-hnu) ) / (2 * hnu);
    if (derivId == id_ncold)
        return dCdnu * nu / ncold[ir];
    real_t lnLEE = lnLambda->evaluateLnLambdaT(ir);
    if (derivId == id_Tcold) {
        real_t dlnLEEdTcold = lnLambda->evaluatePartialAtP(ir, 0, derivId, izs, lnLEE_settings);
        return dCdnu * nu * (dlnLEEdTcold / lnLEE - 2. / Tcold[ir]);
    }

    // derivId == id_ions
    real_t nfree = ions->GetFreeElectronDensityFromQuasiNeutrality(ir);
    if(nfree == 0)
        return 0;
    len_t _, Z0;
    ions->GetIonIndices(izs, _, Z0);
    real_t dZeffdni = Z0/nfree * (Z0 - ions->GetNZ0Z0(ir)/nfree);
    real_t dlnLEEdni = lnLambda->evaluatePartialAtP(ir, 0, derivId, izs, lnLEE_settings);
    real_t hZeff = Zeff * epsilon_forward;
    real_t dCdZeff = ( coefficient(ft[ir], Zeff+hZeff, nu)
                     - coefficient(ft[ir], Zeff, nu) ) / hZeff;

    return dCdnu * nu * (dZeffdni / Zeff + dlnLEEdni / lnLEE) + dCdZeff * dZeffdni;
}

/**
 * Calculates the partial derivative of the coefficient alpha.
 *
 *  ir:       radial cell grid point.
 *  derivId:  unknown quantity ID to differentiate with respect to.
 *  iz:       ion charge state index
 *  iZ:       ion species index
 */
real_t BootstrapCurrent::evaluatePartialCoefficientAlpha(
    len_t ir, len_t derivId, len_t izs, len_t iZ
) {
    if ( (derivId != id_Tcold) && (derivId != id_ions) && (derivId != id_Ni) && (derivId != id_Wi) )
        return 0;
    if ( ((derivId == id_Ni) || (derivId == id_Wi)) && (iZ != iZMain) ) // Only main ion contributes!
        return 0;

    real_t nu = evaluateIonCollisionFrequency(ir);
    real_t Zeff = ions->GetZeff(ir);
    real_t hnu = nu * epsilon_central;
    real_t dAdnu = ( evaluateCoefficientAlpha(ft[ir], Zeff, nu+hnu)
                   - evaluateCoefficientAlpha(ft[ir], Zeff, nu-hnu) ) / (2*hnu);

    real_t lnLII = lnLambda->evaluateLnLambdaII(ir);
    if (derivId == id_Tcold) {
        real_t dlnLIIdTcold = lnLambda->evaluatePartialAtP(ir, 0, derivId, izs, lnLII_settings);
        return dlnLIIdTcold * nu * dAdnu / lnLII;
    }
    // these two only contributes for the main ion, ie. if (iZ == iZMain)
    if (derivId == id_Ni)
        return dAdnu * 3. * nu / NiMain[ir];
    if (derivId == id_Wi)
        return -dAdnu * 2. * Constants::ec * nu / WiMain[ir];

    // derivId == id_ions
    len_t _, Z0;
    ions->GetIonIndices(izs, _, Z0);
    real_t nfree = ions->GetFreeElectronDensityFromQuasiNeutrality(ir);
    real_t dZeffdni = Z0/nfree * (Z0 - ions->GetNZ0Z0(ir)/nfree);
    real_t dlnLIIdni = lnLambda->evaluatePartialAtP(ir, 0, derivId, izs, lnLII_settings);
    real_t hZeff = Zeff * epsilon_forward;
    real_t dAdZeff = ( evaluateCoefficientAlpha(ft[ir], Zeff+hZeff, nu)
                     - evaluateCoefficientAlpha(ft[ir], Zeff, nu) ) / hZeff;

    return dAdnu * nu * (4.* dZeffdni / Zeff + dlnLIIdni / lnLII) + dAdZeff * dZeffdni;
}
