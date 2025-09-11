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
    BootstrapCurrent *bs, IonHandler *ih, len_t iZ, real_t sf
) : BootstrapEquationTerm(g, u, ih, bs, iZ, sf) {

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
real_t BootstrapIonDensityTerm::GetCoefficient(len_t ir) {

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);
    real_t alpha = bs->getCoefficientAlpha(ir);

    real_t coefficient;
    if (bs->includeIonTemperatures)
        coefficient = bs->p[ir] / bs->n[ir] - 2./3. * ( 1. + alpha ) * bs->Wi[rOffset + ir] / (bs->Ni[rOffset + ir] * Constants::ec);
    else
        coefficient = - bs->Tcold[ir] * alpha;

    return pre * l31 * coefficient;
}

/**
 * Partial derivative of Coefficient.
 */
real_t BootstrapIonDensityTerm::GetPartialCoefficient(len_t ir, len_t derivId, len_t jzs, len_t jZ) {

    real_t pre = bs->getConstantPrefactor(ir);
    real_t l31 = bs->getCoefficientL31(ir);
    real_t dl31 = bs->evaluatePartialCoefficientL31(ir, derivId, jzs);
    real_t alpha = bs->getCoefficientAlpha(ir);
    real_t dalpha = bs->evaluatePartialCoefficientAlpha(ir, derivId, jzs, jZ);

    real_t dCoefficient = dl31;
    if (!bs->includeIonTemperatures) {
        dCoefficient *= - bs->Tcold[ir] * alpha;
        dCoefficient -= l31 * bs->Tcold[ir] * dalpha;
        if (derivId == id_Tcold)
            dCoefficient -= l31 * alpha;
    } else {
        dCoefficient *= bs->p[ir] / bs->n[ir] - ( 1. + alpha ) * bs->Wi[rOffset + ir];
        dCoefficient -= l31 * 2./3. * dalpha * bs->Wi[rOffset + ir] / bs->Ni[rOffset + ir];
        if (derivId == id_ncold)
            dCoefficient += l31 * ( bs->Tcold[ir] - bs->p[ir] / bs->n[ir] ) / bs->n[ir];
        else if (derivId == id_Tcold)
            dCoefficient += l31 * bs->ncold[ir] / bs->n[ir];
        else if (derivId == id_Ni)
            dCoefficient += l31 * ( (1. + alpha ) * bs->Wi[rOffset + ir] / ( 1.5 * bs->Ni[rOffset + ir] * bs-> Ni[rOffset + ir]) - bs->p[ir] / (bs->n[ir] * bs->n[ir]) );
        else if (derivId == id_Wi) // IE: Added: l31 / (1.5 * bs->n[ir] * Constants::ec) -
            dCoefficient += l31 / (1.5 * bs->n[ir] * Constants::ec) - l31 * ( 1. + alpha ) / ( 1.5 * bs->Ni[rOffset + ir] * Constants::ec );
    }
    return pre * dCoefficient;
}
