/**
 * Implementation of bootstrap current term driven by the ion density
 * gradients.
 */
#include "DREAM/Equations/Fluid/BootstrapIonDensityTerm.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapIonDensityTerm::BootstrapIonDensityTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u,
    BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, bs, sf) {

    SetUnknownID(id_Ni);

    SetName("BootstrapIonDensityTerm");
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

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);
    real_t alpha = bs->getCoefficientAlpha(ir);

    real_t coefficient = bs->p[ir] / bs->n[ir];
    if (bs->includeIonTemperatures)
        coefficient -= 2./3. * ( 1. + alpha ) * bs->Wi[nr * iZ + ir] / (bs->Ni[nr * iZ + ir] * Constants::ec);

    return pre * l31 * coefficient;
}

/**
 * Partial derivative of Coefficient.
 */
real_t BootstrapIonDensityTerm::GetPartialCoefficient(len_t ir, len_t derivId, len_t index, len_t iZ) {

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);
    real_t dl31 = bs->evaluatePartialCoefficientL31(ir, derivId, index);

    real_t dCoefficient = dl31;
    if (!bs->includeIonTemperatures) {
        dCoefficient *= bs->Tcold[ir];
        if (derivId == id_Tcold)
            dCoefficient += l31;
    } else {
        real_t alpha = bs->getCoefficientAlpha(ir);
        real_t dalpha = bs->evaluatePartialCoefficientAlpha(ir, derivId, index, iZ);

        len_t i = nr * iZ + ir;
        dCoefficient *= bs->p[ir] / bs->n[ir] - ( 1. + alpha ) * bs->Wi[i];
        dCoefficient -= l31 * 2./3. * dalpha * bs->Wi[i] / bs->Ni[i];
        if (derivId == id_ncold)
            dCoefficient += l31 * ( bs->Tcold[ir] - bs->p[ir] / bs->n[ir] ) / bs->n[ir];
        else if (derivId == id_Tcold)
            dCoefficient += l31 * bs->ncold[ir] / bs->n[ir];
        else if (derivId == id_Ni)
            dCoefficient += l31 * ( (1. + alpha ) * bs->Wi[i] / ( 1.5 * bs->Ni[i] * bs-> Ni[i]) - bs->p[ir] / (bs->n[ir] * bs->n[ir]) );
        else if (derivId == id_Wi)
            dCoefficient += l31 * ( 1. + alpha ) / ( 1.5 * bs->Ni[i] * Constants::ec );
    }
    return pre * dCoefficient;
}
