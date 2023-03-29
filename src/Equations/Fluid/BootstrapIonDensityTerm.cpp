/**
 * Implementation of bootstrap current term driven by the ion density
 * gradients.
 */
#include "DREAM/Equations/Fluid/BootstrapIonDensityTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapIonDensityTerm::BootstrapIonDensityTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, bs, sf) {

    SetUnknownID(id_Ni);

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
real_t BootstrapIonDensityTerm::GetCoefficient(len_t ir, len_t iZ) {
    real_t coefficient = bs->p[ir] / bs->n[ir];
    if (bs->includeIonTemperatures)
        coefficient -= 2./3. * (1. + bs->evaluateCoefficientAlpha(ir)) * bs->Wi[nr * iZ + ir] / bs->Ni[nr * iZ + ir];
    return bs->constantPrefactor[ir] * bs->coefficientL31[ir] * coefficient;
}

/**
 * Partial derivative of Coefficient.
 */
real_t BootstrapIonDensityTerm::GetPartialCoefficient(
    len_t ir, len_t derivId, len_t index, len_t iZ
) {
    real_t dCoefficient = bs->evaluatePartialCoefficientL31(ir, derivId, index);
    if (!bs->includeIonTemperatures) {
        dCoefficient *= bs->Tcold[ir];
        dCoefficient += bs->coefficientL31[ir] * (real_t)(derivId == id_Tcold);
    } else {
        dCoefficient *= bs->p[ir] / bs->n[ir];
        real_t dAlpha = bs->evaluatePartialCoefficientAlpha(ir, derivId, index, iZ);
        len_t i = nr * iZ + ir;
        dCoefficient -= bs->coefficientL31[ir] * 2./3. * dAlpha * bs->Wi[i] / bs->Ni[i];
        if (derivId == id_ncold)
            dCoefficient += bs->coefficientL31[ir] * (bs->Tcold[ir] - bs->p[ir] / bs->n[ir]) / bs->n[ir];
        else if (derivId == id_Tcold)
            dCoefficient += bs->coefficientL31[ir] * bs->ncold[ir] / bs->n[ir];
        else if (derivId == id_Ni) {
            real_t alpha = bs->evaluateCoefficientAlpha(ir);
            dCoefficient += bs->coefficientL31[ir] * ( 1. / (bs->n[ir] * bs->n[ir]) + (1. + alpha) * 2./3. * bs->Wi[i] / (bs->Ni[i] *bs-> Ni[ir]) );
        } else if (derivId == id_Wi) {
            real_t alpha = bs->evaluateCoefficientAlpha(ir);
            dCoefficient += bs->coefficientL31[ir] * 2./3. * ( 1./ bs->n[ir] - (1. + alpha) / bs->Ni[i] );
        }
    }
    return bs->constantPrefactor[ir] * dCoefficient;
}
