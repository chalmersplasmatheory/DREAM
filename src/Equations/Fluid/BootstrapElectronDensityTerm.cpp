/**
 * Implementation of bootstrap current term driven by the electron density
 * gradient.
 */
#include "DREAM/Equations/Fluid/BootstrapElectronDensityTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapElectronDensityTerm::BootstrapElectronDensityTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, BootstrapCurrent *bs, IonHandler *ih, real_t sf
) : BootstrapEquationTerm(g, u, ih, sf), Bootstrap(bs) {

    id_X = id_ncold;

    SetName("BootstrapElectronDensityTerm");
    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Tcold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_ions);
    if (includeIonTemperatures)
        AddUnknownForJacobian(u, id_Wi);
}

/**
 * Coefficient = A * L31 * p / n
 */
real_t BootstrapElectronDensityTerm::GetCoefficient(len_t ir, len_t /* iZ */) {
    real_t coefficient;
    if (includeIonTemperatures)
        coefficient = p[ir] / n[ir];
    else
        coefficient = Tcold[ir];
    return constantPrefactor[ir] * coefficientL31[ir] * coefficient;
}

real_t BootstrapElectronDensityTerm::GetPartialCoefficient(
    len_t ir, len_t derivId, len_t index, len_t /* iZ */
) {
    real_t dCoefficient = evaluatePartialCoefficientL31(ir, derivId, index)
    if (!includeIonTemperatures)
        dCoefficient *= (real_t)(derivId == id_Tcold);
    else {
        dCoefficient *= p[ir] / n[ir];
        switch (derivId) {
            case id_ncold:
                dCoefficient += coefficientL31[ir] * (Tcold[ir] - p[ir] / n[ir]) / n[ir];
                break;
            case id_Tcold:
                dCoefficient += coefficientL31[ir] * ncold[ir] / n[ir];
                break;
            case id_Ni:
                dCoefficient -= coefficientL31[ir] / (n[ir] * ni[ir]);
                break;
            case id_Wi:
                dCoefficient += coefficientL31[ir] * 2. / (3. * n[ir]);
        }
    }
    return constantPrefactor[ir] * dCoefficient;
}
