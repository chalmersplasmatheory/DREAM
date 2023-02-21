/**
 * Implementation of bootstrap current term driven by the ion thermal energy
 * gradients.
 */
#include "DREAM/Equations/Fluid/BootstrapIonThermalEnergyTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapIonThermalEnergyTerm::BootstrapIonThermalEnergyTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, sf), Bootstrap(bs) {

    id_X = id_Wi;
    iZMainIon = iZMain;

    SetName("BootstrapIonThermalEnergyTerm");
    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Tcold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_ions);
}

/**
 * If ion temperature evolution is included, then for each ion:
 *      Coefficient = A * L31 * (1 + alpha) * 2/3
 *
 * If not included, ie. Ti = Tcold, then this term is not used!
 */
real_t BootstrapIonThermalEnergyTerm::GetCoefficient(len_t ir, len_t iZ) {
    real_t coefficient = p[ir] / n[ir] * 2./3. * (1. + evaluateCoefficientAlpha(ir));
    return constantPrefactor[ir] * coefficientL31[ir] * coefficient;
}

/**
 * Partial derivative of Coefficient.
 */
real_t BootstrapIonThermalEnergyTerm::GetPartialCoefficient(
    len_t ir, len_t derivId, len_t index, len_t iZ
) {
    real_t dCoefficient = evaluatePartialCoefficientL31(ir, derivId, index);
    dCoefficient *= 1. + evaluateCoefficientAlpha(ir);
    dCoefficient += coefficientL31[ir] * evaluatePartialCoefficientAlpha(ir, derivId, index, iZ);
    return constantPrefactor[ir] * 2./3. * dCoefficient;
}
