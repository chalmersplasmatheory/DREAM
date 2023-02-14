/**
 * Implementation of a class that calculates and stores quantities related to
 * bootstrap current. For further details on the implementation, see doc/notes/bootstrap.
 */
 #include "DREAM/Equations/BootstrapCurrent.hpp"

 /**
  * Constructor.
  */
BootstrapCurrent::BootstrapCurrent(
    FVM::RadialGrid *rg, FMV::UnknownQuantityHandler *u,
    IonHandler *ih, CoulombLogarithm *lnLee, CoulombLogarithm *lnLii,
) {

    rGrid = rg;
    unknowns = u;
    ions = ih;
    lnLambdaEE = lnLee;
    lnLambdaii = lnLii;


    id_jtot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);

    // Memory allocation
    AllocateQuantities();

    // small number for numerical differentiation
    epsilon = sqrt(std::numeric_limits<real_t>::epsilon());

    // equilibrium constants
    const real_t R0 = grid->GetR0();
    for (len_t ir = 0; ir < nr; ir++) {
        // calculate the geometric prefactor
        const real_t *BtorGOverR0 = rGrid->GetBTorG(ir);        // G / R0
        const real_t *FSA_B2 = rGrid->GetFSA_B2(ir);            // <B^2> / Bmin^2
        const real_t *Bmin = rGrid->GetBmin(ir);                // Bmin
        const real_t *psiPrimeRef = rGrid->GetPsiPrimeRef(ir);  // R0 d(psi_ref)/dr
        constantPrefactor[ir] = - BtorGOverR0 * R0 * R0 / ( FSA_B2 * Bmin * psiPrimeRef);
        if (ir == 0 || ir == nr-1) // boundary cases
            constantPrefactor[ir] /= 2 * rGrid->GetDr( (ir != 0)*(nr-2) );
        else
            constantPrefactor[ir] /= ( rGrid->GetDr(ir-1) + rGrid->GetDr(ir+1) );

        // calculate geometric constants
        ft[ir] = 1 - rGrid->GetEffPassFrac(ir);
    }

    // locate the main ion index
    for (len_t iZ = 0; iZ < ions.GetNZ(); iZ++)
        if (ions.GetZ(iZ) == 1) {
            iZMain = iZ;
            break;
        }
}

 /**
  * Destructor.
  */
BootstrapCurrent::~BootstrapCurrent() {
    DeallocateQuantities();
}

/**
 * Allocate memory for arrays stored in this object.
 */
void BootstrapCurrent::AllocateQuantities() {
    DeallocateQuantities();

    constantPrefactor   = new real_t[nr];
    coefficientL31      = new real_t[nr];
    Zeff                = new real_t[nr];
    ft                  = new real_t[nr];
    qR0                 = new real_t[nr];
    n                   = new real_t[nr];
    p                   = new real_t[nr];
}

/**
 * Deallocate memory for arrays stored in this object.
 */
void BootstrapCurrent::DeallocateQuantities() {
    delete [] constantPrefactor;
    delete [] coefficientL31;
    delete [] Zeff;
    delete [] ft;
    delete [] qR0;
    delete [] n;
    delete [] p;
}



void BootstrapCurrent::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) {
    // currentX = u->GetUnknownData(id_X);
    // previousX = u-GetUnknownDataPrevious(id_X);

    *jtot  = unknowns->GetUnknownData(id_jtot);
    *ncold = unknowns->GetUnknownData(id_ncold);
    *ni    = unknowns->GetUnknownData(id_ni);
    *Tcold = unknowns->GetUnknownData(id_Tcold);
    *Wi    = unknowns->GetUnknownData(id_Wi);




    for (len_t ir = 0; ir < nr; ir++) {

        // calculate safety factor (normalised to R0), used in the collision frequencies
        real_t Ip = TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir, rGrid, jtot);
        qR0[ir] = fabs(grid->SafetyFactorNormalized(ir, Constants::mu0 * Ip));

        // calculate effective charge
        Zeff[ir] = ions->GetZeff(ir);

        // calculate L31 as it is used in all terms
        coefficientL31[ir] = evaluateCoefficientL31(ir);

        // obtain main ion data (relevant only for alpha coefficient)
        niMain[ir] = ions->GetTotalIonDensity(ir, iZMain);
        WiMain[ir] = Wi[nr * ions->GetIndex(iZMain, 1) + ir];

        // calculate total pressure and number density
        n = ncold[ir];
        p = ncold[ir] * Tcold[ir];
        for (iZ = 0; iZ < ions->GetNZ(); iZ++) {
            n += ions->GetTotalIonDensity(ir, iZ);
            p += 2. * Wi[nr * ions->GetIndex(iZ, 1) + ir] / 3.;
        }
    }
}


/**
 * Thermal electron collision frequency (as defined in Eq. 18b(d) in Sauter et al. 1999).
 */
real_t BootstrapCurrent::evaluateElectronCollisionFrequency(len_t ir, real_t ncold, real_t Tcold, real_t Zeff) {
    real_t lnLee = lnLambdaEE->evaluateLnLambdaT();
    real_t eps = rGrid->GetR(ir) / rGrid->GetR0();
    return 6.921e-18 * ncold * lnLee * Zeff * qR0[ir] / (eps * sqrt(eps) * Tcold * Tcold);
}

real_t BootstrapCurrent::evaluateElectronCollisionFrequency(len_t ir) {
    return evaluateElectronCollisionFrequency(ir, ncold[ir], Tcold[ir], Zeff[ir]);
}

/**
 * Ion collision frequency (SORT OF as defined in Eq. 18c(e) in Sauter et al. 1999).
 */
real_t BootstrapCurrent::evaluateIonCollisionFrequency(len_t ir, real_t ni, real_t Wi, real_t Zeff) {
    real_t lnLii = lnLambdaII->evaluateLnLambdaII(ir); // formula from Wesson
    real_t eps = rGrid->GetR(ir) / rGrid->GetR0();
    return 4.90e-18 * 9. * ni*ni*ni * lnLii * Zeff*Zeff*Zeff*Zeff * qR0[ir] / (eps * sqrt(eps) * 4. * Wi*Wi );
}


/**
 * Main ion collision frequency, used for all ion species in Redl's model
 */
real_t BootstrapCurrent::evaluateIonCollisionFrequency(len_t ir) {
    return evaluateIonCollisionFrequency(ir, niMain[ir], WiMain[ir], Zeff[ir]);
}


// CALCULATIONS OF THE BOOTSTRAP COEFFICIENTS


/**
 * Calculates the coefficient L31 (as defined in Eqs. 10-11 in Redl et al. 2021).
 */
real_t BootstrapCurrent::evaluateCoefficientL31(len_t ir, real_t ncold, real_t Tcold, real_t Zeff, real_t nu) {

    // Eq. 11
    real_t d = (.52 + .086 * sqrt(nu)) * (1. + .87 * ft[ir]) * nu / (1. + 1.13 * sqrt(Zeff - 1.));
    d += 1 + .67*(1. - .7*ft) * sqrt(nu) / (.56 + .44*Zeff);
    real_t f31 = ft[ir] / d;

    // Eq. 10
    real_t a = 1 / (pow(Zeff, 1.2) - .71);
    return (1. + .15*a)*f31 - .22*a*f31*f31 + .01*a*f31*f31*f31 + .06*a*f31*f31*f31*f31;
}

real_t BootstrapCurrent::evaluateCoefficientL31(len_t ir) {
    real_t nu = evaluateElectronCollisionFrequency(ir, ncold[ir], Tcold[ir], Zeff[ir]);
    return evaluateCoefficientL31(ir, ncold[ir], Tcold[ir], Zeff[ir], nu);
}

/**
 * Calculates the coefficient L32 (as defined in Eqs. 12-16 in Redl et al. 2021).
 */
real_t BootstrapCurrent::evaluateCoefficientL32(len_t ir, real_t ncold, real_t Tcold, real_t Zeff, real_t nu) {

    // Eq. 14
    real_t dee = sqrt( 1. + 2. * sqrt(Zeff - 1.) );
    dee += ft[ir] * ft[ir] * sqrt( nu * (.075 + .25 * (Zeff - 1.) * ( Zeff - 1.)) );
    dee *= .13 * (1. - .38 * ft[ir]) * nu / (Zeff * Zeff);
    dee += 1. + .23 * (1. - .96 * ft[ir]) * sqrt(nu / Zeff);
    real_t f32ee = ft[ir] / dee;

    // Eq. 13
    real_t F32ee = (.1 + .6 * Zeff) * f32ee * (1. - f32ee * f32ee * f32ee);
    F32ee /= Zeff * (.77 + .63*(1. + pow(Zeff-1., 1.1)));
    F32ee += .7 / (1. + .2 * Zeff) * f32ee * f32ee * (1. - 1.3 * f32ee + .2 * f32ee * f32ee);
    F32ee += 1.3 * f32ee * f32ee * f32ee * f32ee / (1. + .5 * Zeff);

    // Eq. 16
    real_t dei = 1. + .87 * (1. + .39 * ft[ir]) * sqrt(nu) / (1. + 2.95 * (Zeff - 1.) * (Zeff - 1.));
    dei += 1.53 * (1. - .37 * ft[ir]) * nu * (2. + .375 * (Zeff - 1.));
    real_t f32ei = ft[ir] / dei;

    // Eq. 15
    real_t F32ei = 5.5 * f32ei * f32ei * (1. - .8 * f32ei - .2 * f32ei * f32ei);
    F32ei /= (1.5 + 2 * Zeff);
    F32ei -= (.4 + 1.93 * Zeff) * f32ei * (1. - f32ei * f32ei * f32ei) / (Zeff * (.8 + .6 * Zeff));
    F32ei -= 1.3 * f32ei * f32ei * f32ei * f32ei / (1. + .5 * Zeff);

    return F32ee + F32ei;
}

real_t BootstrapCurrent::evaluateCoefficientL32(len_t ir) {
    real_t nu = evaluateElectronCollisionFrequency(ir, ncold[ir], Tcold[ir], Zeff[ir]);
    return evaluateCoefficientL32(ir, ncold[ir], Tcold[ir], Zeff[ir], nu);
}

/**
 * Calculates the coefficient alpha (as defined in Eqs. 20-21 in Redl et al. 2021).
 */
real_t BootstrapCurrent::evaluateCoefficientAlpha(len_t ir, real_t Zeff, real_t nu) {

    // Eq. 20
    real_t alpha0 = - ( .6 + .055 * (Zeff - 1.) ) / ( .53 + .17 * (Zeff - 1.) );
    real_t ft2 = ft[ir] * ft[ir];
    alpha0 *= (1. - ft[ir]) / (1. - ( .31 - .065 * (Zeff - 1.)) ft[ir] - .25 * ft2);

    // Eq. 21
    real_t alpha = ( alpha0 + .7 * Zeff sqrt(ft[ir] * nu) ) / ( 1. + .18 * sqrt(nu) );
    real_t ft6 = ft2 * ft2 * ft2;
    alpha -= .002 * nu*nu * ft6;
    alpha /= 1. + .004 * nu*nu * ft6;
    return alpha;
}

// NEEDS TO BE REVISED
real_t BootstrapCurrent::evaluateCoefficientAlpha(len_t ir) {
    real_t nu = evaluateIonCollisionFrequency(ir);
    return evaluateCoefficientAlpha(ir, Zeff[ir], nu);
}


// CALCULATIONS OF PARTIAL DERIVATIVES OF THE BOOTSTRAP COEFFICIENTS

/**
 * Calculates the partial derivate of the coefficient L31.
 */
real_t BootstrapCurrent::evaluatePartialCoefficientL31(len_t ir, len_t derivId, len_t index) {
    return evaluateNumericalDerivative(ir, derivId, index, &evaluateCoefficientL32);
}


/**
 * Calculates the partial derivate of the coefficient L32.
 */
real_t BootstrapCurrent::evaluatePartialCoefficientL32(len_t ir, len_t derivId, len_t index) {
    return evaluateNumericalDerivative(ir, derivId, index, &evaluateCoefficientL32);
}


/**
 * Calculates the partial derivative of coefficient L31 or L32.
 */
real_t BootstrapCurrent::evaluateNumericalDerivative(
    len_t ir, len_t derivId, len_t index,
real_t (*coefficient)(len_t, real_t, real_t, real_t, real_t)
) {
    real_t nu = evaluateElectronCollisionFrequency(ir);
    if (derivId == id_ncold) {
        real_t h = ncold[ir] * epsilon;
        return ( coefficient(ir, ncold[ir]+h, Tcold[ir], Zeff[ir], nu)
               - coefficient(ir, ncold[ir]-h, Tcold[ir], Zeff[ir], nu) ) / (2*h);
    } else if (derivId == id_Tcold) {
        real_t h = Tcold[ir] * epsilon;
        return ( coefficient(ir, ncold[ir], Tcold[ir]+h, Zeff[ir], nu)
               - coefficient(ir, ncold[ir], Tcold[ir]-h, Zeff[ir], nu) ) / (2*h);
    } else if (derivId == id_ni) {
        len_t iZ, Z0;
        ions->GetIonIndices(index, iZ, Z0);
        real_t nfree = ions->GetFreeElectronDensityFromQuasiNeutrality(ir);
        // if(nfree==0) // why would we need this?
        //     return 0;
        real_t dZeffdni = Z0/nfree * (Z0 - ions->GetNZ0Z0(ir)/nfree);
        real_t dnudZeff = nu / Zeff[ir];
        real_t hnu = nu * epsilon;
        real_t dL31dnu = ( coefficient(ir, ncold[ir], Tcold[ir], Zeff[ir], nu+hnu)
                         - coefficient(ir, ncold[ir], Tcold[ir], Zeff[ir], nu-hnu) ) / (2*hnu));
        real_t hZeff = Zeff[ir] * epsilon;
        real_t dL31dZeff = ( coefficient(ir, ncold[ir], Tcold[ir], Zeff[ir]+hZeff, nu)
                           - coefficient(ir, ncold[ir], Tcold[ir], Zeff[ir]-hZeff, nu) ) / (2*hZeff);
        return dZeffdni * (dnudZeff * dL31dnu + dL31dZeff);
    } else
        return 0;
}

 /**
  * Calculates the partial derivate of the coefficient alpha. (YET TO BE IMPLEMENTED!!)
  */
real_t BootstrapCurrent::evaluatePartialCoefficientAlpha(len_t ir, len_t derivId, len_t index) {
    real_t nu = evaluateIonCollisionFrequency(ir);
    } if (derivId == id_Tcold) {
        real_t dlnLiidTcold = 1. / Tcold[ir];
        real_t dnudlnLii = nu / lnLambdaII->evaluateLnLambdaII(ir);
        real_t hnu = nu * epsilon;
        real_t dalphadnu = ( evaluateCoefficientAlpha(ir, Zeff[ir], nu-hnu)
                           - evaluateCoefficientAlpha(ir, Zeff[ir], nu+hnu) ) / (2*hnu);
        return dlnLiidTcold * dnudlnLii * dalphadnu;
    } else if (derivId == id_ni) {
        len_t iZ, Z0;
        ions->GetIonIndices(index, iZ, Z0);
        real_t nfree = ions->GetFreeElectronDensityFromQuasiNeutrality(ir);
        real_t dZeffdni = Z0/nfree * (Z0 - ions->GetNZ0Z0(ir)/nfree);
        real_t dnudZeff = 4. * nu / Zeff[ir];
        real_t hnu = nu * epsilon;
        real_t dalphadnu = ( evaluateCoefficientAlpha(ir, Zeff[ir], nu-hnu)
                           - evaluateCoefficientAlpha(ir, Zeff[ir], nu+hnu) ) / (2*hnu);
        real_t hZeff = Zeff[ir] * epsilon;
        real_t dalphadZeff = ( evaluateCoefficientAlpha(ir, Zeff[ir]-hZeff, nu)
                             - evaluateCoefficientAlpha(ir, Zeff[ir]+hZeff, nu) ) / (2*hZeff);
        real_t dlnLiidni = - .5 * Z0 / nfree;
        real_t dnudlnLii = nu / lnLambdaII->evaluateLnLambdaII(ir);
        real_t dnudni = dZeffdni * dnudZeff + dlnLiidni * dnudlnLii;
        if (iZ == iZMain)
            dnudni += nu * ( 3. / niMain[ir] - 2. / WiMain[ir] );
        return dZeffdni * dalphadZeff + dnudni * dalphadnu;
    } else
        return 0;
}
