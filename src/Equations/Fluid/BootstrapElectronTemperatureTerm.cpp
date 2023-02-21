/**
 * Implementation of bootstrap current term driven by the electron temperature
 * gradient.
 */
#include "DREAM/Equations/Fluid/BootstrapElectronTemperatureTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapElectronTemperatureTerm::BootstrapElectronTemperatureTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, sf), Bootstrap(bs) {

    id_X = id_Tcold;

    SetName("BootstrapElectronTemperatureTerm");
    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Tcold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_ions);
}

/**
 * If ion temperature evolution is included:
 *      Coefficient = A * (L31 + L32) * ncold
 *
 * If not included, ie. Ti = Tcold, then:
 *      Coefficient = A * [ (L31 + L32) * ncold + L31 * (1 + alpha) * ionsum(ni) ]
 */
real_t BootstrapElectronTemperatureTerm::GetCoefficient(len_t ir, len_t /* iZ */) {
    real_t coefficient = ( coefficientL31[ir] + evaluateCoefficientL32(ir) ) * ncold[ir];
    if (!includeIonTemperatures) {
        real_t nitot = 0;
        for (i = ir; i < nr * nZ; i += nr)
            nitot += Ni[i];
        coefficient += coefficientL31[ir] * (1. + evaluateCoefficientAlpha(ir)) * nitot;
    }
    return constantPrefactor[ir] * coefficient;
}

/**
 * Partial derivative of Coefficient.
 */
real_t BootstrapElectronTemperatureTerm::GetPartialCoefficient(
    len_t ir, len_t derivId, len_t index, len_t iZ
) {
    real_t dl31 = evaluatePartialCoefficientL31(ir, derivId, index);
    real_t dl32 = evaluatePartialCoefficientL32(ir, derivId, index);
    real_t dCoefficient = (dl31 + dl32) * ncold[ir];
    if (derivId == id_ncold)
        dCoefficient += coefficientL31[ir] + evaluateCoefficientL32(ir);
    if (!includeIonTemperatures) {
        real_t nitot = 0;
        for (i = ir; i < nr * nZ; i += nr)
            nitot += Ni[i];
        real_t alpha = evaluateCoefficientAlpha(ir);
        dCoefficient += dl31 * (1. + alpha) * nitot;

        real_t dAlpha = evaluatePartialCoefficientAlpha(ir, derivId, index, iZ);
        dCoefficient += oefficientL31[ir] * dAlpha * nitot;

        if (derivId == id_Ni)
            dCoefficient += coefficientL31[ir] * (1. + alpha);
    }
    return constantPrefactor[ir] * dCoefficient;
}
