/**
 * Implementation of bootstrap current term driven by the ion density
 * gradients.
 */
#include "DREAM/Equations/Fluid/BootstrapElectronDensityTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapElectronDensityTerm::BootstrapElectronDensityTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, sf), Bootstrap(bs) {

    id_X = id_Ni;

    SetName("BootstrapElectronDensityTerm");
    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Tcold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_ions);
}

/**
 * If ion temperature evolution is included, then for each ion:
 *      Coefficient = A * L31 * [ p / n - (1 + alpha) * 2/3 * Wi / ni ]
 *
 * If not included, ie. Ti = Tcold, then for each ion:
 *      Coefficient = A * L31 *  p / n
 */
real_t BootstrapElectronDensityTerm::GetCoefficient(len_t ir, len_t iZ) {
    real_t coefficient = p[ir] / n[ir];
    if (includeIonTemperatures)
        coefficient -= 2./3. * (1. + evaluateCoefficientAlpha(ir)) * Wi[nr * iZ + ir] / Ni[nr * iZ + ir] );
    return constantPrefactor[ir] * coefficientL31[ir] * coefficient;
}

/**
 * Partial derivative of Coefficient.
 */
real_t BootstrapElectronDensityTerm::GetPartialCoefficient(
    len_t ir, len_t derivId, len_t index, len_t iZ
) {
    real_t dCoefficient = evaluatePartialCoefficientL31(ir, derivId, index);
    if (!includeIonTemperatures)
        dCoefficient *= (real_t)(derivId == id_Tcold);
    else {
        dCoefficient *= p[ir] / n[ir];
        real_t dAlpha = evaluatePartialCoefficientAlpha(ir, derivId, index, iZ);
        len_t i = nr * iZ + ir;
        dCoefficient -= coefficientL31[ir] * 2./3. * dAlpha * Wi[i] / Ni[i];
        switch (derivId) {
            case id_ncold:
                dCoefficient += coefficientL31[ir] * (Tcold[ir] - p[ir] / n[ir]) / n[ir];
                break;
            case id_Tcold:
                dCoefficient += coefficientL31[ir] * ncold[ir] / n[ir];
            case id_Ni:
                real_t alpha = evaluateCoefficientAlpha(ir);
                dCoefficient += coefficient31[ir] * ( 1. / (n[ir] * n[ir]) + (1. + alpha) * 2./3. * Wi[i] / (Ni[i] * Ni[ir]) );
                break;
            case id_Wi:
                real_t alpha = evaluateCoefficientAlpha(ir);
                dCoefficient += coefficient31[ir] * 2./3. * ( 1./ n[ir] - (1. + alpha) / Ni[i] );
                break;
        }
    }
    return constantPrefactor[ir] * dCoefficient;
}
